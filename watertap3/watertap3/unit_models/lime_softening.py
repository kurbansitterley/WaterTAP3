from copy import deepcopy
from pyomo.environ import Var, Param, check_optimal_termination, units as pyunits
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

import idaes.logger as idaeslog
from idaes.core import declare_process_block_class
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.constants import Constants
from idaes.core.util.scaling import set_scaling_factor, get_scaling_factor

from watertap.costing.util import make_capital_cost_var
from watertap3.core.wt3_unit_sido import WT3UnitProcessSIDOData


## REFERENCE
## CAPITAL:
# citation here
## ELECTRICITY:
# citation here

module_name = "lime_softening"
# Costing model adapted from McGivney & Kawamura (2007)
# Figure 2.4.1b + Appendix A5a
# only including: mix > lime addition > flocculator (solids contact) > CO2 contact tank > sedimentation/clarifer
rapid_mix_capital_coeffs = {
    300: (1268.1, 0.4431),
    600: (1131.2, 0.4731),
    900: (618.01, 0.5696),
}


def cost_lime_softening(blk):
    blk.basis_year = 2007
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")

    blk.rapid_mix_capital_cost_base = Var(
        initialize=rapid_mix_capital_coeffs[blk.unit_model.rapid_mix_velocity_grad][0],
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Rapid mix capital cost base",
    )

    blk.rapid_mix_capital_cost_exp = Var(
        initialize=rapid_mix_capital_coeffs[blk.unit_model.rapid_mix_velocity_grad][1],
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Rapid mix capital cost exponent",
    )

    blk.lime_feed_system_capital_cost_base = Var(
        initialize=12985,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Lime feed system capital cost base",
    )

    blk.lime_feed_system_capital_cost_exp = Var(
        initialize=0.5901,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Lime feed system capital cost exponent",
    )

    blk.floc_20_capital_cost_base = Var(
        initialize=908675,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Floc 20 capital cost base",
    )

    blk.floc_20_capital_cost_exp = Var(
        initialize=0.7191,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Floc 20 capital cost exponent",
    )

    blk.floc_50_capital_cost_base = Var(
        initialize=1000977,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Floc 50 capital cost base",
    )

    blk.floc_50_capital_cost_exp = Var(
        initialize=0.7579,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Floc 50 capital cost exponent",
    )

    blk.floc_80_capital_cost_base = Var(
        initialize=1218085,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Floc 80 capital cost base",
    )

    blk.floc_80_capital_cost_exp = Var(
        initialize=0.8266,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Floc 80 capital cost exponent",
    )

    blk.CO2_reactor_capital_cost_base = Var(
        initialize=3470.6,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="CO2 reactor capital cost base",
    )

    blk.CO2_reactor_capital_cost_exp = Var(
        initialize=0.6173,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="CO2 reactor capital cost exp",
    )

    blk.sedimentation_basin_capital_cost_base = Var(
        initialize=13572,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Sedimentation basin capital cost base",
    )

    blk.sedimentation_basin_capital_cost_exp = Var(
        initialize=0.3182,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Sedimentation basin capital cost exp",
    )

    blk.number_redundant_train = Var(
        initialize=1,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Number of redundant treatment trains",
    )

    blk.pump_efficiency = Var(
        initialize=0.8,
        bounds=(0, 1),
        units=pyunits.dimensionless,
        doc="Pump efficiency",
    )
    blk.motor_efficiency = Var(
        initialize=0.8,
        bounds=(0, 1),
        units=pyunits.dimensionless,
        doc="Motor efficiency",
    )

    blk.lift_height = Var(
        initialize=100,
        bounds=(1, 1e5),
        units=pyunits.feet,
        doc="Pump lift height for chemical flows",
    )

    blk.handle_costing_unit_params()
    blk.fix_all_vars()

    blk.rapid_mix_capital_cost = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Rapid mix capital cost",
    )

    blk.flocculator_capital_cost = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Flocculator capital cost",
    )

    blk.CO2_reactor_capital_cost = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="CO2 reactor capital cost",
    )

    blk.sedimentation_basin_capital_cost = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Sedimentation basin capital cost",
    )

    blk.lime_feed_system_capital_cost = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Lime feed system capital cost",
    )
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TPEC")

    capital_cost_expr = 0

    @blk.Constraint(doc="Rapid mix capital cost")
    def rapid_mix_capital_cost_constraint(b):
        rapid_mix_vol_dim = pyunits.convert(
            b.unit_model.rapid_mix_volume / pyunits.gallon,
            to_units=pyunits.dimensionless,
        )
        return (
            b.rapid_mix_capital_cost
            == b.rapid_mix_capital_cost_base
            * rapid_mix_vol_dim**b.rapid_mix_capital_cost_exp
        )

    capital_cost_expr += blk.unit_model.number_rapid_mix * blk.rapid_mix_capital_cost

    @blk.Constraint(
        doc="Flocculator capital cost"
    )  # assumes 3 stage flocculator: 20 > 50 > 80 [vel_gradient]
    def flocculator_capital_cost_constraint(b):
        floc_vol_dim = pyunits.convert(
            b.unit_model.flocculator_volume / pyunits.Mgallon,
            to_units=pyunits.dimensionless,
        )
        return (
            b.flocculator_capital_cost
            == b.floc_20_capital_cost_base * floc_vol_dim**b.floc_20_capital_cost_exp
            + b.floc_50_capital_cost_base * floc_vol_dim**b.floc_50_capital_cost_exp
            + b.floc_80_capital_cost_base * floc_vol_dim**b.floc_80_capital_cost_exp
        )

    capital_cost_expr += (
        blk.unit_model.number_flocculator * blk.flocculator_capital_cost
    )

    @blk.Constraint(doc="CO2 reactor capital cost")
    def CO2_reactor_capital_cost_constraint(b):
        CO2_reactor_area_dim = pyunits.convert(
            b.unit_model.CO2_reactor_area / pyunits.ft**2,
            to_units=pyunits.dimensionless,
        )
        return (
            b.CO2_reactor_capital_cost
            == b.CO2_reactor_capital_cost_base
            * CO2_reactor_area_dim**b.CO2_reactor_capital_cost_exp
        )

    capital_cost_expr += (
        blk.unit_model.number_CO2_reactor * blk.CO2_reactor_capital_cost
    )

    @blk.Constraint(doc="Sedimentation basin capital cost")
    def sedimentation_basin_capital_cost_constraint(b):
        sed_area_dim = pyunits.convert(
            b.unit_model.sedimentation_surface_area / pyunits.ft**2,
            to_units=pyunits.dimensionless,
        )
        return (
            b.sedimentation_basin_capital_cost
            == b.sedimentation_basin_capital_cost_base
            * sed_area_dim**b.sedimentation_basin_capital_cost_exp
        )

    capital_cost_expr += (
        blk.unit_model.number_sedimentation * blk.sedimentation_basin_capital_cost
    )

    @blk.Constraint(doc="Lime feed system capital cost")
    def lime_feed_system_capital_cost_constraint(b):
        mass_flow_lime_dim = pyunits.convert(
            b.unit_model.mass_flow_lime / (pyunits.lb / pyunits.day),
            to_units=pyunits.dimensionless,
        )
        return (
            b.lime_feed_system_capital_cost
            == b.lime_feed_system_capital_cost_base
            * mass_flow_lime_dim**b.lime_feed_system_capital_cost_exp
        )

    capital_cost_expr += (
        blk.unit_model.number_lime_feed * blk.lime_feed_system_capital_cost
    )

    @blk.Expression(doc="Capital cost of redundant train")
    def capital_cost_redundant_train(b):
        return b.number_redundant_train * (
            b.rapid_mix_capital_cost
            + b.flocculator_capital_cost
            + b.CO2_reactor_capital_cost
            + b.sedimentation_basin_capital_cost
            + b.lime_feed_system_capital_cost
        )

    capital_cost_expr += blk.capital_cost_redundant_train

    @blk.Constraint(doc="Capital cost constraint")
    def capital_cost_constraint(b):
        return b.capital_cost == pyunits.convert(
            capital_cost_expr, to_units=b.costing_package.base_currency
        )

    pumping_power_expr = 0

    @blk.Expression(doc="Pumping power for lime suspension")
    def pumping_power_lime_suspension(b):
        return pyunits.convert(
            (
                blk.unit_model.mass_flow_lime
                * Constants.acceleration_gravity
                * blk.lift_height
            )
            / (blk.pump_efficiency * blk.motor_efficiency),
            to_units=pyunits.kW,
        )

    pumping_power_expr += blk.pumping_power_lime_suspension

    @blk.Expression(doc="Pumping power for lime suspension")
    def pumping_power_soda_ash_suspension(b):
        return pyunits.convert(
            (
                blk.unit_model.mass_flow_soda_ash
                * Constants.acceleration_gravity
                * blk.lift_height
            )
            / (blk.pump_efficiency * blk.motor_efficiency),
            to_units=pyunits.kW,
        )

    pumping_power_expr += blk.pumping_power_soda_ash_suspension

    @blk.Expression(doc="Floc paddle power")
    def pumping_power_floc(b):
        floc_power_expr = 0
        floc_vol = pyunits.convert(
            b.unit_model.flocculator_volume, to_units=pyunits.m**3
        )
        for g in [20, 50, 80]:
            g = g * pyunits.s**-1
            floc_power_expr += (
                pyunits.convert(
                    (g**2 * floc_vol * b.unit_model.properties_in.visc_d)
                    / b.motor_efficiency,
                    to_units=pyunits.kilowatt,
                )
                * b.unit_model.number_flocculator
            )
        return floc_power_expr

    pumping_power_expr += blk.pumping_power_floc

    @blk.Expression(doc="Mixer power")
    def pumping_power_mixer(b):
        rapid_mix_vol = pyunits.convert(
            b.unit_model.rapid_mix_volume, to_units=pyunits.m**3
        )
        g = b.unit_model.rapid_mix_velocity_grad * pyunits.s**-1
        return (
            pyunits.convert(
                (g**2 * rapid_mix_vol * b.unit_model.properties_in.visc_d)
                / b.motor_efficiency,
                to_units=pyunits.kilowatt,
            )
            * b.unit_model.number_rapid_mix
        )

    pumping_power_expr += blk.pumping_power_mixer

    blk.cost_chemical_flow(material="lime_suspension", dose=blk.unit_model.lime_dose)
    blk.cost_chemical_flow(material="soda_ash", dose=blk.unit_model.soda_ash_dose)
    blk.costing_package.cost_flow(pumping_power_expr, "electricity")


@declare_process_block_class("LimeSoftening")
class UnitProcessData(WT3UnitProcessSIDOData):
    def build(self):
        super().build()

        if "rapid_mix_velocity_grad" in self.config.unit_params.keys():
            self.rapid_mix_velocity_grad = self.config.unit_params[
                "rapid_mix_velocity_grad"
            ]
        else:
            self.rapid_mix_velocity_grad = 900

        self.rapid_mix_retention_time = Param(
            initialize=10,
            mutable=True,
            units=pyunits.second,
            doc="Rapid mix retention time",
        )

        self.floc_retention_time = Param(
            initialize=10,
            mutable=True,
            units=pyunits.minute,
            doc="Flocculator retention time",
        )

        self.CO2_reactor_retention_time = Param(
            initialize=45,
            mutable=True,
            units=pyunits.second,
            doc="CO2 reactor retention time",
        )

        self.CO2_reactor_height = Param(
            initialize=10,
            mutable=True,
            units=pyunits.ft,
            doc="CO2 reactor height",
        )

        self.sedimentation_loading_rate = Param(
            initialize=0.75,  # between 0.5-1 gpm/ft2
            mutable=True,
            units=pyunits.gallon / pyunits.minute / pyunits.ft**2,
            doc="Loading rate for sedimentation basins",
        )

        self.lime_dose = Param(
            initialize=0.005,
            mutable=True,
            units=pyunits.kg / pyunits.m**3,
            doc="Lime dose",
        )

        self.lime_solution_density = Param(
            initialize=1250,
            mutable=True,
            units=pyunits.kg / pyunits.m**3,
            doc="Lime solution density",
        )

        self.soda_ash_dose = Param(
            initialize=0,
            mutable=True,
            units=pyunits.kg / pyunits.m**3,
            doc="Soda ash dose",
        )

        self.soda_ash_solution_density = Param(
            initialize=1250,
            mutable=True,
            units=pyunits.kg / pyunits.m**3,
            doc="Soda ash solution density",
        )

        self.rapid_mix_volume = Var(
            initialize=500,
            bounds=(0, 145000),
            units=pyunits.gallon,
            doc="Rapid mix basin volume",
        )

        self.number_rapid_mix = Var(
            initialize=1,
            bounds=(1, None),
            units=pyunits.dimensionless,
            doc="Number rapid mixers",
        )

        self.flocculator_volume = Var(
            initialize=0.5,
            bounds=(0, None),
            units=pyunits.Mgallon,
            doc="Flocculator basin volume",
        )

        self.number_flocculator = Var(
            initialize=2,
            bounds=(1, None),
            units=pyunits.dimensionless,
            doc="Number flocculator",
        )

        self.CO2_reactor_area = Var(
            initialize=1000,
            bounds=(0, 32300),
            units=pyunits.ft**2,
            doc="CO2 reactor floor area",
        )

        self.number_CO2_reactor = Var(
            initialize=1,
            bounds=(1, None),
            units=pyunits.dimensionless,
            doc="Number CO2 reactors",
        )

        self.number_lime_feed = Var(
            initialize=2,
            bounds=(1, None),
            units=pyunits.dimensionless,
            doc="Number lime feed system",
        )

        self.sedimentation_surface_area = Var(
            initialize=10000,
            units=pyunits.ft**2,
            bounds=(0, 150000),
            doc="Basin surface area",
        )

        self.number_sedimentation = Var(
            initialize=1,
            bounds=(1, None),
            units=pyunits.dimensionless,
            doc="Number sedimentation basins",
        )

        self.handle_unit_params()

        @self.Expression(doc="Mass flow lime")
        def mass_flow_lime(b):
            return pyunits.convert(
                b.lime_dose * b.properties_in.flow_vol,
                to_units=pyunits.lb / pyunits.day,
            )

        @self.Expression(doc="Mass flow soda ash")
        def mass_flow_soda_ash(b):
            return pyunits.convert(
                b.soda_ash_dose * b.properties_in.flow_vol,
                to_units=pyunits.lb / pyunits.day,
            )

        @self.Expression(doc="Total rapid mix volume calculation")
        def total_rapid_mix_volume(b):
            return pyunits.convert(
                b.properties_in.flow_vol * b.rapid_mix_retention_time,
                to_units=pyunits.gallon,
            )

        @self.Constraint(doc="Rapid mix volume per reactor")
        def rapid_mix_volume_constraint(b):
            return b.rapid_mix_volume == b.total_rapid_mix_volume / b.number_rapid_mix

        @self.Expression(doc="Total flocculator volume required")
        def total_flocculator_volume(b):
            return pyunits.convert(
                b.properties_in.flow_vol * b.floc_retention_time,
                to_units=pyunits.Mgallon,
            )

        @self.Constraint(doc="Flocculator volume per basin")
        def floc_volume_constraint(b):
            return (
                b.flocculator_volume
                == b.total_flocculator_volume / b.number_flocculator
            )

        @self.Expression(doc="Total CO2 reactor volume required")
        def total_CO2_reactor_volume(b):
            return pyunits.convert(
                b.properties_in.flow_vol * b.CO2_reactor_retention_time,
                to_units=pyunits.gallon,
            )

        @self.Expression(doc="Total CO2 reactor area required")
        def total_CO2_reactor_area(b):
            return pyunits.convert(
                b.total_CO2_reactor_volume / b.CO2_reactor_height,
                to_units=pyunits.ft**2,
            )

        @self.Constraint(doc="Total CO2 reactor area")
        def CO2_reactor_area_constraint(b):
            return b.CO2_reactor_area == b.total_CO2_reactor_area / b.number_CO2_reactor

        @self.Expression(doc="Total sedimentation basin area required")
        def total_sedimentation_basin_area(b):
            return pyunits.convert(
                b.properties_in.flow_vol / b.sedimentation_loading_rate,
                to_units=pyunits.ft**2,
            )

        @self.Constraint(doc="Toal sedimentation basin surface area")
        def sedimentation_surface_area_constraint(b):
            return (
                b.sedimentation_surface_area
                == b.total_sedimentation_basin_area / b.number_sedimentation
            )

        @self.Constraint(doc="Number of rapid mix is half the number of flocculators")
        def number_rapid_mix_floc(b):
            return b.number_rapid_mix == b.number_flocculator

        @self.Constraint(doc="Number of flocculators is equal to number CO2 reactors")
        def number_floc_CO2(b):
            return b.number_flocculator == b.number_CO2_reactor

        @self.Constraint(
            doc="Number of lime feed systems is equal to number flocculators"
        )
        def number_floc_lime(b):
            return b.number_flocculator == b.number_lime_feed

        @self.Constraint(
            doc="Number sedimentation basins is twice the number flocculators"
        )
        def number_floc_sed(b):
            return b.number_sedimentation == b.number_flocculator

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
                state_args_out[k] = v * self.water_recovery.value
            elif k == "conc_mass_comp":
                state_args_out[k] = dict()
                for j, u in v.items():
                    state_args_out[k][j] = (1 - self.removal_fraction[j].value) * u
            elif k == "flow_mass_comp":
                state_args_out[k] = dict()
                for j, u in v.items():
                    if j == "H2O":
                        state_args_out[k][j] = self.water_recovery.value * u
                    else:
                        state_args_out[k][j] = (1 - self.removal_fraction[j].value) * u

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        self.state_args_waste = state_args_waste = deepcopy(state_args)

        for k, v in state_args.items():
            if k == "flow_vol":
                state_args_waste[k] = v * self.water_recovery.value
            elif k == "conc_mass_comp":
                state_args_waste[k] = dict()
                for j, u in v.items():
                    state_args_waste[k][j] = self.removal_fraction[j].value * u
            elif k == "flow_mass_comp":
                state_args_waste[k] = dict()
                for j, u in v.items():
                    if j == "H2O":
                        state_args_waste[k][j] = (1 - self.water_recovery.value) * u
                    else:
                        state_args_waste[k][j] = self.removal_fraction[j].value * u

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

    @property
    def default_costing_method(self):
        return cost_lime_softening
