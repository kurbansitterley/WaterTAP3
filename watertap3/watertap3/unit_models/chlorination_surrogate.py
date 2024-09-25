import os
import sys
from copy import deepcopy
from io import StringIO

from pyomo.environ import (
    Param,
    Var,
    PositiveReals,
    check_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

import idaes.logger as idaeslog
from idaes.core import declare_process_block_class
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.exceptions import InitializationError
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
from idaes.core.util.scaling import set_scaling_factor, get_scaling_factor

from watertap.costing.util import make_capital_cost_var

from watertap3.core.wt3_unit_siso import WT3UnitProcessSISOData

## REFERENCE:
# CAPITAL: Table 3.23 - User's Manual for Integrated Treatment Train Toolbox - Potable Reuse (IT3PR) Version 2.0

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
__author__ = "Kurban Sitterley"

surrogate_file = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/surrogates/chlorination_surrogate.json"
)

module_name = "chlorination"


def cost_chlorination(blk):

    blk.basis_year = 2014
    blk.basis_currency = getattr(pyunits, f"kUSD_{blk.basis_year}")

    blk.energy_intensity = Var(
        initialize=5e-5,
        bounds=(0, None),
        units=pyunits.kWh / pyunits.m**3,
        doc="Chlorination energy intensity",
    )

    blk.fix_all_vars()

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)

    blk.flow_mgd = Var(
        initialize=1,
        bounds=(0, None),
        units=pyunits.Mgallons / pyunits.day,
        doc="Flow in MGD for surrogate model",
    )

    blk.capital_cost_surrogate = Var(
        initialize=1e3,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Capital cost from surrogate model",
    )

    stream = StringIO()
    oldstdout = sys.stdout
    sys.stdout = stream

    blk.surrogate_blk = SurrogateBlock(concrete=True)
    blk.surrogate = PysmoSurrogate.load_from_file(surrogate_file)
    blk.surrogate_blk.build_model(
        blk.surrogate,
        input_vars=[blk.unit_model.dose, blk.flow_mgd],
        output_vars=[blk.capital_cost_surrogate],
    )

    sys.stdout = oldstdout

    @blk.Constraint(doc="Flow in MGD conversion for surrogate model")
    def flow_mgd_constraint(b):
        return b.flow_mgd == pyunits.convert(
            b.unit_model.properties_in.flow_vol,
            to_units=pyunits.Mgallons / pyunits.day,
        )
    
    @blk.Constraint(doc="Capital cost equation")
    def capital_cost_constraint(b):
        return b.capital_cost == pyunits.convert(
            b.capital_cost_surrogate, to_units=b.costing_package.base_currency
        )

    @blk.Expression(doc="Power required")
    def power_required(b):
        return pyunits.convert(
            b.energy_intensity * b.unit_model.properties_in.flow_vol,
            to_units=pyunits.kilowatt,
        )

    blk.costing_package.cost_flow(blk.power_required, "electricity")
    blk.cost_chemical_flow()


@declare_process_block_class("Chlorination")
class UnitProcessData(WT3UnitProcessSISOData):
    def build(self):
        super().build()

        if "chemical" in self.config.unit_params.keys():
            chemical = self.config.unit_params["chemical"]
        else:
            chemical = "chlorine"
        self.properties_in.flow_vol

        self.contact_time = Param(
            initialize=1.5,
            mutable=True,
            domain=PositiveReals,
            units=pyunits.hour,
            doc="Contact time",
        )
        self.ct = Param(
            initialize=450,
            mutable=True,
            domain=PositiveReals,
            units=(pyunits.mg * pyunits.minute) / pyunits.liter,
            doc="CT value",
        )
        self.decay_rate = Param(
            initialize=5,
            mutable=True,
            domain=PositiveReals,
            units=pyunits.mg / (pyunits.liter * pyunits.hour),
            doc="Chlorine decay rate",
        )

        self.dose = Var(
            initialize=5,
            bounds=(0.1, 26),
            units=pyunits.mg / pyunits.liter,
            doc="Chlorine dose",
        )

        self.handle_unit_params()
        self.chemical = chemical

        if not self.dose.is_fixed():

            @self.Constraint(doc="Chlorine dose equation")
            def eq_chlorine_dose(b):
                dose_a = b.decay_rate * b.contact_time
                dose_b = pyunits.convert(
                    b.ct / b.contact_time, to_units=pyunits.mg / pyunits.liter
                )
                return b.dose == pyunits.convert(
                    dose_a + dose_b, to_units=pyunits.mg / pyunits.liter
                )

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for initialization routines

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        # return
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        if solver is None:
            opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        flags = self.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info("Initialization Step 1a Complete.")

        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # Set state_args from inlet state
        if state_args is None:
            self.state_args = state_args = {}
            state_dict = self.properties_in.define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        self.state_args_out = state_args_out = deepcopy(state_args)
        for k, v in state_args.items():
            if k == "flow_vol":
                state_args_out[k] = v
            elif k == "conc_mass_comp":
                state_args_out[k] = dict()
                for j, u in v.items():
                    state_args_out[k][j] = (1 - self.removal_fraction[j].value) * u

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        if not self.dose.is_fixed():
            cvc(self.dose, self.eq_chlorine_dose)
            
        if hasattr(self, "costing"):
            cvc(self.costing.flow_mgd, self.costing.flow_mgd_constraint)
            self.costing.flow_mgd.fix()
            self.costing.flow_mgd_constraint.deactivate()
            cvc(
                self.costing.capital_cost_surrogate,
                self.costing.surrogate_blk.pysmo_constraint["capital_cost"],
            )
            self.costing.capital_cost_surrogate.fix()
            self.costing.surrogate_blk.pysmo_constraint["capital_cost"].deactivate()
            cvc(self.costing.capital_cost, self.costing.capital_cost_constraint)

            init_log.info("Initialization Step 1c Complete.")

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)
                res = opt.solve(self, tee=slc.tee)

        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Release Inlet state
        self.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize.")
        
        self.initialized = True


    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        set_scaling_factor(self.dose, value(self.dose)**-1)

    @property
    def default_costing_method(self):
        return cost_chlorination