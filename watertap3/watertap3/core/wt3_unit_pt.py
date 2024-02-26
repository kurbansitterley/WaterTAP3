# “Institute for the Design of Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE Framework) Copyright (c) 2019, by the software owners:
# The Regents of the University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University Research Corporation, et al.
# All rights reserved."
import idaes.logger as idaeslog
from idaes.core import declare_process_block_class
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.scaling import set_scaling_factor, get_scaling_factor
from pyomo.environ import check_optimal_termination, Var, value, units as pyunits
from pyomo.network import Port
from .wt3_unit_base import WT3UnitProcessBaseData

__all__ = ["WT3UnitProcessPT"]

__author__ = "Kurban Sitterley"

solver = get_solver()


@declare_process_block_class("WT3UnitProcessPT")
class WT3UnitProcessPTData(WT3UnitProcessBaseData):
    """
    WaterTAP3 passthrough unit.
    Unit models using this class have no component removal and no waste stream.
    """

    def build(self):
        super().build()
        units_meta = self.config.property_package.get_metadata().get_derived_units

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True
        self.properties_in = prop_in = self.config.property_package.state_block_class(
            doc="Material properties of inlet stream", **tmp_dict
        )
        tmp_dict["defined_state"] = False
        self.properties_out = prop_out = self.config.property_package.state_block_class(
            doc="Material properties of outlet stream", **tmp_dict
        )

        @self.Constraint(
            self.config.property_package.component_list, doc="Component mass balances"
        )
        def component_mass_balance(b, j):
            return prop_in.flow_mass_comp[j] == prop_out.flow_mass_comp[j]

        self.inlet = Port(noruleinit=True, doc="Inlet Port")
        self.inlet.add(prop_in.flow_mass_comp, "flow_mass")
        # self.inlet.add(prop_in.flow_vol, 'flow_vol')
        # self.inlet.add(prop_in.conc_mass_comp, 'conc_mass')
        # self.inlet.add(prop_in.temperature, 'temperature')
        # self.inlet.add(prop_in.pressure, 'pressure')

        self.outlet = Port(noruleinit=True, doc="Outlet Port")
        self.outlet.add(prop_out.flow_mass_comp, "flow_mass")
        # self.outlet.add(prop_out.flow_vol, 'flow_vol')
        # self.outlet.add(prop_out.conc_mass_comp, 'conc_mass')
        # self.outlet.add(prop_out.temperature, 'temperature')
        # self.outlet.add(prop_out.pressure, 'pressure')

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

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info("Initialization Step 1b Complete.")

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
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
