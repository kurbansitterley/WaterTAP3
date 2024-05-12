from copy import deepcopy
from pyomo.environ import (
    Param,
    NonNegativeReals,
    check_optimal_termination,
    Var,
    value,
    units as pyunits,
)
from pyomo.network import Port
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

import idaes.logger as idaeslog
from idaes.core import declare_process_block_class
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.scaling import set_scaling_factor, get_scaling_factor

from watertap.costing.util import make_capital_cost_var, make_fixed_operating_cost_var

from watertap3.core.wt3_unit_sido import WT3UnitProcessSIDOData
from watertap.core.util.model_diagnostics.infeasible import *

_log = idaeslog.getLogger(__name__)

module_name = "reverse_osmosis"

__author__ = "Kurban Sitterley"


def cost_reverse_osmosis(blk):

    blk.basis_year = 2007
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")

    blk.number_trains = Var(
        initialize=2,
        units=pyunits.dimensionless,
        doc="Membrane cost",
    )
    blk.membrane_unit_cost = Var(
        initialize=30,
        bounds=(10, 80),
        units=blk.basis_currency / pyunits.m**2,
        doc="Membrane cost",
    )
    blk.vessel_unit_cost = Var(
        initialize=1000,
        units=blk.basis_currency,
        doc="Unit cost for pressure vessel",
    )
    blk.rack_unit_cost = Var(
        initialize=33,
        units=blk.basis_currency / pyunits.ft,
        doc="Unit cost for rack per length",
    )
    blk.number_vessels_per_area = Var(
        initialize=0.025,
        units=pyunits.m**-2,
        doc="Number of pressure vessels per unit membrane area",
    )
    blk.rack_length_base = Var(
        initialize=150,
        units=pyunits.ft,
        doc="Base length for racks",
    )
    blk.rack_length_additional = Var(
        initialize=5,
        units=pyunits.ft,
        doc="Length required per additional vessel",
    )
    blk.factor_membrane_replacement = Var(
        initialize=0.25,
        bounds=(0.01, 1),
        units=pyunits.year**-1,
        doc="Membrane replacement rate",
    )
    blk.factor_chemical_cost = Var(
        initialize=0.01,
        bounds=(0.001, 0.05),
        units=pyunits.year**-1,
        doc="Reverse osmosis chemical addition costs",
    )

    blk.pump_capital_cost_base = Var(
        initialize=1.908,  # 53 / 1e5 * 3600,
        bounds=(1, 3),
        units=blk.basis_currency,  # ($ * hr) / (m^3 * bar)
        doc="Pump capital cost equation - base",
    )
    blk.pump_capital_cost_exp = Var(
        initialize=0.97,
        bounds=(0, 1),
        units=pyunits.dimensionless,
        doc="Pump capital cost equation - exponent",
    )
    blk.pump_efficiency = Var(
        initialize=0.8,
        bounds=(0, 1),
        units=pyunits.dimensionless,
        doc="Pump efficiency",
    )

    blk.erd_capital_cost_base = Var(
        initialize=3134.7,
        units=blk.basis_currency,
        doc="ERD capital cost equation - base",
    )
    blk.erd_capital_cost_exp = Var(
        initialize=0.58,
        bounds=(0, 1),
        units=pyunits.dimensionless,
        doc="ERD capital cost equation - exponent",
    )
    blk.erd_efficiency = Var(
        initialize=0.8,
        bounds=(0, 1),
        units=pyunits.dimensionless,
        doc="Energy recovery device efficiency",
    )

    blk.fix_all_vars()

    make_fixed_operating_cost_var(blk)
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    blk.pressure_vessel_capital_cost = Var(
        initialize=1e6,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Pressure vessel cost",
    )
    blk.rack_support_capital_cost = Var(
        initialize=1e6,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Rack support cost",
    )
    blk.pump_capital_cost = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Pump capital cost",
    )
    blk.erd_capital_cost = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Energy recovery device capital cost",
    )
    blk.membrane_capital_cost = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Membrane capital cost",
    )

    @blk.Expression(doc="Number of pressure vessels per train")
    def number_vessels(b):
        return blk.unit_model.membrane_area * b.number_vessels_per_area

    @blk.Expression(doc="Total length for rack support per train")
    def rack_support_length(b):
        return b.rack_length_base + b.number_vessels * b.rack_length_additional

    @blk.Expression(doc="Total capital cost of pressure vessels and rack.")
    def support_capital_cost(b):
        return 3.3 * (b.pressure_vessel_capital_cost + b.rack_support_capital_cost)

    @blk.Expression(doc="Pump power")
    def pump_power(b):
        return pyunits.convert(
            (
                (b.unit_model.operating_pressure - b.unit_model.pressure_atm)
                * b.unit_model.properties_in.flow_vol
            )
            / b.pump_efficiency,
            to_units=pyunits.watt,
        )

    if blk.unit_model.has_erd:

        @blk.Expression(doc="ERD power")
        def erd_power(b):
            return pyunits.convert(
                (
                    (b.unit_model.properties_waste.pressure - b.unit_model.pressure_atm)
                    * b.unit_model.properties_waste.flow_vol
                )
                / b.erd_efficiency,
                to_units=pyunits.watt,
            )

        @blk.Constraint(doc="ERD capital cost equation")
        def erd_capital_cost_constraint(b):
            flow_mass_waste = pyunits.convert(
                b.unit_model.properties_waste.flow_mass_comp["tds"]
                + b.unit_model.properties_waste.flow_mass_comp["H2O"],
                to_units=pyunits.kg / pyunits.hr,
            )
            flow_mass_waste_dim = pyunits.convert(
                flow_mass_waste * pyunits.hr / pyunits.kg,
                to_units=pyunits.dimensionless,
            )
            return (
                b.erd_capital_cost
                == b.erd_capital_cost_base
                * flow_mass_waste_dim**b.erd_capital_cost_exp
            )

    else:

        @blk.Expression(doc="ERD power")
        def erd_power(b):
            return 0 * pyunits.watt

        blk.erd_capital_cost.fix(0)

    @blk.Expression(doc="Total power required")
    def power_required(b):
        return pyunits.convert(
            b.pump_power, to_units=pyunits.kilowatt
        ) - pyunits.convert(b.erd_power, to_units=pyunits.kilowatt)

    @blk.Constraint(doc="Membrane capital cost equation")
    def membrane_capital_cost_constraint(b):
        return (
            b.membrane_capital_cost == b.unit_model.membrane_area * b.membrane_unit_cost
        )

    @blk.Constraint(doc="Pump capital cost equation")
    def pump_capital_cost_constraint(b):
        pump_power_dim = pyunits.convert(
            b.pump_power / pyunits.watt, to_units=pyunits.dimensionless
        )
        return (
            b.pump_capital_cost
            == b.pump_capital_cost_base * pump_power_dim**b.pump_capital_cost_exp
        )

    @blk.Constraint(doc="Pressure vessel capital cost equation")
    def pressure_vessel_capital_cost_constraint(b):
        return (
            b.pressure_vessel_capital_cost
            == b.number_vessels * b.vessel_unit_cost * b.number_trains
        )

    @blk.Constraint(doc="Rack capital cost equation")
    def rack_support_capital_cost_constraint(b):
        return (
            b.rack_support_capital_cost
            == b.rack_support_length * b.rack_unit_cost * b.number_trains
        )

    @blk.Constraint(doc="Capital cost equation")
    def capital_cost_constraint(b):
        return b.capital_cost == pyunits.convert(
            b.pump_capital_cost
            + b.membrane_capital_cost
            + b.erd_capital_cost
            + b.support_capital_cost,
            to_units=b.costing_package.base_currency,
        )

    @blk.Constraint(doc="Fixed operating cost equation")
    def fixed_operating_cost_constraint(b):
        return b.fixed_operating_cost == pyunits.convert(
            (b.factor_chemical_cost * b.capital_cost),
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        ) + pyunits.convert(
            (
                b.factor_membrane_replacement
                * b.membrane_unit_cost
                * b.unit_model.membrane_area
            ),  # assume chemical addition costs are 1% of capital cost
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )

    blk.costing_package.cost_flow(blk.power_required, "electricity")
    # TODO: add chemical costing = 1% of capital costs
    # self.chem_dict = {"unit_cost": 0.01}


@declare_process_block_class("ReverseOsmosis")
class UnitProcessData(WT3UnitProcessSIDOData):
    def build(self):

        super().build()

        self.properties_in.pressure.set_value(3e6)
        self.properties_waste.pressure.set_value(5e6)

        if (
            "erd" not in self.config.unit_params.keys()
            or not self.config.unit_params["erd"]
        ):
            self.has_erd = False
        else:
            self.has_erd = True

        self.pressure_atm = Param(
            initialize=101325,
            mutable=True,
            units=pyunits.Pa,
            doc="Atmospheric pressure",
        )

        self.deltaP_waste = Param(
            initialize=3e5,
            units=pyunits.Pa,
            doc="Pressure change between inlet and waste",
        )

        self.module_membrane_area = Param(
            initialize=37.16,
            mutable=True,
            units=pyunits.m**3,
            doc="Membrane area per module",
        )

        self.water_permeability = Param(
            initialize=4.2e-12,
            mutable=True,
            units=pyunits.m / (pyunits.Pa * pyunits.s),
            doc="Water permeability",
        )

        self.salt_permeability = Param(
            initialize=3.5e-8,
            mutable=True,
            units=pyunits.m / pyunits.s,
            doc="Salt permeability",
        )

        self.target_water_recovery = Param(
            initialize=0.5,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Target water recovery fraction",
        )

        self.target_permeate_salinity = Param(
            initialize=0.5,
            mutable=True,
            units=pyunits.kg / pyunits.m**3,
            doc="Target permeate salinity",
        )

        self.target_tds_removal_fraction = Param(
            initialize=0.98,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Target TDS removal fraction",
        )

        self.deltaP_outlet = Var(
            initialize=1e-6,
            units=pyunits.Pa,
            doc="Pressure change between inlet and outlet",
        )

        self.membrane_area = Var(
            initialize=1e6,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Total membrane area",
        )

        self.operating_pressure = Var(
            initialize=1e6,
            domain=NonNegativeReals,
            bounds=(101324, 8.3e6),
            units=pyunits.Pa,
            doc="Operating pressure",
        )

        @self.Expression(doc="Number membrane modules needed")
        def number_modules(b):
            return b.membrane_area / b.module_membrane_area

        # @self.Expression(doc="")
        # def operating_pressure(b):
        #     return (b.properties_in.pressure + b.properties_waste.pressure) * 0.5

        @self.Expression(doc="")
        def avg_osmotic_pressure(b):
            return (
                b.properties_in.pressure_osmotic + b.properties_waste.pressure_osmotic
            ) * 0.5

        @self.Constraint(doc="Operating presssure")
        def eq_operating_pressure(b):
            return (
                b.operating_pressure
                == (b.properties_in.pressure + b.properties_waste.pressure) * 0.5
            )

        # @self.Constraint(doc="Permeate pressure")
        # def eq_permeate_pressure(b):
        #     return b.deltaP_outlet == -b.properties_in.pressure + b.pressure_atm

        @self.Constraint(doc="Water transport equation")
        def eq_water_transport(b):
            return (
                b.properties_out.flow_vol
                == b.water_permeability
                # * (b.properties_in.pressure - b.properties_in.pressure_osmotic)
                * (b.operating_pressure - b.avg_osmotic_pressure) * b.membrane_area
            )

        @self.Constraint(doc="Salt transport equation")
        def eq_salt_transport(b):
            return (
                b.properties_out.flow_mass_comp["tds"]
                == b.salt_permeability
                * (
                    b.properties_in.conc_mass_comp["tds"]
                    + b.properties_waste.conc_mass_comp["tds"]
                    # - b.properties_out.conc_mass_comp["tds"]
                )
                * b.membrane_area
            )

        @self.Constraint()
        def eq_tds_removal_lb(b):
            return b.removal_fraction["tds"] >= b.target_tds_removal_fraction

        @self.Constraint()
        def eq_pressure_lb(b):
            return b.operating_pressure >= 1.01 * b.avg_osmotic_pressure

        @self.Expression()
        def operating_pressure_bar(b):
            return pyunits.convert(b.operating_pressure, to_units=pyunits.bar)

        @self.Expression()
        def avg_osmotic_pressure_bar(b):
            return pyunits.convert(b.avg_osmotic_pressure, to_units=pyunits.bar)

        # @self.Constraint(doc="Pressure of brine stream")
        # def eq_pressure_waste(b):
        #     return (
        #         b.properties_waste.pressure == b.properties_in.pressure + b.deltaP_waste
        #     )

        # @self.Constraint(doc="Pressure of perm stream")
        # def eq_pressure_permeate(b):
        #     return (
        #         b.properties_out.pressure == b.properties_in.pressure + b.deltaP_outlet
        #     )

        # @self.Constraint(doc="Inlet pressure must be greater than osmotic")
        # def eq_min_pressure_inlet(b):
        #     return b.properties_in.pressure >= 1.05* b.properties_in.pressure_osmotic

        # @self.Constraint(doc="Outlet pressure must be greater than osmotic")
        # def eq_min_pressure_brine(b):
        #     return b.properties_waste.pressure >= 1.05* b.properties_waste.pressure_osmotic

        # @self.Constraint(doc="Target water recovery")
        # def eq_target_water_recovery_lb(b):
        #     return b.water_recovery >= 0.9 * b.target_water_recovery

        # @self.Constraint(doc="Target water recovery")
        # def eq_target_water_recovery_ub(b):
        #     return b.water_recovery <= 1.05 * b.target_water_recovery

        # @self.Constraint(doc="Target permeate salinity")
        # def eq_target_perm_salinity(b):
        #     return b.properties_out.conc_mass_comp["tds"] <= b.target_permeate_salinity

        self.handle_unit_params()

    def initialize(self, **kwargs):
        self.initialize_build(**kwargs)

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
        self.removal_fraction["tds"].unfix()

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
            if k == "flow_mass_comp":
                for j, u in v.items():
                    if j == "H2O":
                        state_args_out[k][j] = u * self.target_water_recovery.value
                    else:
                        rf = self.flowsheet().unit_removal_fractions[self.unit_name][j]
                        state_args_out[k][j] = u * (1 - rf)
            #     state_args_out[k] == v * 0.5
            # elif k == "conc_mass_comp":
            #     state_args_out[k] == dict()
            #     for j, u in v.items():
            #         state_args_out[k][j] = (1 - self.removal_fraction[j].value) * u

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        self.state_args_waste = state_args_waste = deepcopy(state_args)
        for k, v in state_args.items():
            if k == "flow_mass_comp":
                for j, u in v.items():
                    if j == "H2O":
                        state_args_waste[k][j] = u * (
                            1 - self.target_water_recovery.value
                        )
                    else:
                        rf = self.flowsheet().unit_removal_fractions[self.unit_name][j]
                        state_args_waste[k][j] = u * rf

        self.properties_waste.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_waste,
        )

        init_log.info("Initialization Step 1c Complete.")

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

        print_infeasible_constraints(self)
        # self.removal_fraction["tds"].unfix()

        # if not check_optimal_termination(res):
        #     raise InitializationError(f"Unit model {self.name} failed to initialize.")
        self.initialized = True

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for p in [self.properties_in, self.properties_out, self.properties_waste]:
            set_scaling_factor(p.pressure, 1e-5)
            set_scaling_factor(p.pressure_osmotic, 1e-5)
        set_scaling_factor(self.operating_pressure, 1e-5)
        set_scaling_factor(self.membrane_area, 1e-4)
        set_scaling_factor(self.deltaP_outlet, 1e-5)

    @property
    def default_costing_method(self):
        return cost_reverse_osmosis
