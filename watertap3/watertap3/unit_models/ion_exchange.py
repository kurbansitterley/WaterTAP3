import math
from copy import deepcopy

from pyomo.environ import (
    Var,
    value,
    Param,
    check_optimal_termination,
    units as pyunits,
)
import idaes.logger as idaeslog

from idaes.core import declare_process_block_class
from idaes.core.util.constants import Constants
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.exceptions import InitializationError, ConfigurationError

from watertap.costing.util import make_capital_cost_var, make_fixed_operating_cost_var
from watertap3.core.wt3_unit_siso import WT3UnitProcessSISOData
from watertap3.core.wt3_unit_sido import WT3UnitProcessSIDOData
from watertap3.core.util.pumping_energy import pumping_energy

## REFERENCE
## CAPITAL:
# Parameter data for all costing variables taken from EPA-WBS cation + anion exchange models
# https://www.epa.gov/sdwa/drinking-water-treatment-technology-unit-cost-models
## ELECTRICITY:
# pumping energy for influent flow

module_name = "ion_exchange"


def cost_ion_exchange(blk):
    blk.basis_year = 2020
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")

    col_vol_gal = pyunits.convert(blk.unit_model.column_volume, to_units=pyunits.gallon)
    bed_vol_ft3 = pyunits.convert(blk.unit_model.bed_volume, to_units=pyunits.ft**3)
    regen_tank_vol_gal = pyunits.convert(
        blk.unit_model.regeneration_volume, to_units=pyunits.gallon
    )
    rinse_backwashing_vol_gal = pyunits.convert(
        blk.unit_model.rinse_volume + blk.unit_model.backwashing_volume,
        to_units=pyunits.gallon,
    )

    blk.anion_exchange_resin_cost = Var(
        initialize=256,
        units=blk.basis_currency / pyunits.ft**3,
        doc="Anion exchange resin cost per cubic ft. Assumes strong base polystyrenic macroporous Type II. From EPA-WBS cost model.",
    )
    blk.cation_exchange_resin_cost = Var(
        initialize=231,
        units=blk.basis_currency / pyunits.ft**3,
        doc="Cation exchange resin cost per cubic ft. Assumes strong acid polystyrenic macroporous. From EPA-WBS cost model.",
    )
    # Ion exchange pressure vessels costed with power equation, col_vol in gallons:
    #   pressure_vessel_cost = A * col_vol ** b

    blk.vessel_A_coeff = Var(
        initialize=1596.499333,
        units=blk.basis_currency,
        doc="Ion exchange pressure vessel cost equation - A coeff., Carbon steel w/ stainless steel internals",
    )
    blk.vessel_b_coeff = Var(
        initialize=0.459496809,
        units=pyunits.dimensionless,
        doc="Ion exchange pressure vessel cost equation - b coeff., Carbon steel w/ stainless steel internals",
    )

    # Ion exchange backwash/rinse tank costed with power equation, tank_vol in gallons:
    #   bw_tank_cost = A * tank_vol ** b

    blk.backwash_tank_A_coeff = Var(
        initialize=308.9371309,
        units=blk.basis_currency,
        doc="Ion exchange backwash tank cost equation - A coeff., Steel tank",
    )
    blk.backwash_tank_b_coeff = Var(
        initialize=0.501467571,
        units=pyunits.dimensionless,
        doc="Ion exchange backwash tank cost equation - b coeff., Steel tank",
    )
    # Ion exchange regeneration solution tank costed with power equation, tank_vol in gallons:
    #   regen_tank_cost = A * tank_vol ** b

    blk.regen_tank_A_coeff = Var(
        initialize=57.02158923,
        units=blk.basis_currency,
        doc="Ion exchange regen tank cost equation - A coeff. Stainless steel",
    )
    blk.regen_tank_b_coeff = Var(
        initialize=0.729325391,
        units=pyunits.dimensionless,
        doc="Ion exchange regen tank cost equation - b coeff. Stainless steel",
    )
    blk.annual_resin_replacement_factor = Var(
        initialize=0.05,
        bounds=(0, None),
        units=pyunits.year**-1,
        doc="Fraction of ion excange resin replaced per year, 4-5% of bed volume - EPA",
    )

    blk.handle_costing_unit_params()
    blk.add_pumping_energy()
    blk.fix_all_vars()
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    make_fixed_operating_cost_var(blk)

    blk.capital_cost_anion_exchange_resin = Var(
        initialize=1e4,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Capital cost for anion exchange resin for one column",
    )
    blk.capital_cost_cation_exchange_resin = Var(
        initialize=1e4,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Capital cost for cation exchange resin for one column",
    )
    blk.capital_cost_vessel = Var(
        initialize=1e4,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Capital cost for one vessel (column)",
    )
    blk.capital_cost_regen_tank = Var(
        initialize=1e4,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Capital cost for regeneration tank",
    )
    blk.capital_cost_backwash_tank = Var(
        initialize=1e4,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Capital cost for backwashing + rinse tank",
    )
    blk.resin_replacement_operating_cost = Var(
        initialize=1e4,
        bounds=(0, None),
        units=blk.basis_currency / blk.costing_package.base_period,
        doc="Operating cost for resin replacement",
    )

    @blk.Constraint(doc="Capital cost for anion exchange resin - single column")
    def capital_cost_anion_exchange_resin_constraint(b):
        return (
            b.capital_cost_anion_exchange_resin
            == bed_vol_ft3 * b.anion_exchange_resin_cost
        )

    @blk.Constraint(doc="Capital cost for anion exchange resin - single column")
    def capital_cost_cation_exchange_resin_constraint(b):
        return (
            b.capital_cost_cation_exchange_resin
            == bed_vol_ft3 * b.cation_exchange_resin_cost
        )

    @blk.Constraint(doc="Capital cost for single vessel")
    def capital_cost_vessel_constraint(b):
        col_vol_gal_dim = pyunits.convert(
            col_vol_gal * pyunits.gal**-1, to_units=pyunits.dimensionless
        )
        return (
            b.capital_cost_vessel
            == b.vessel_A_coeff * col_vol_gal_dim**b.vessel_b_coeff
        )

    @blk.Constraint(doc="Capital cost of regeneration tank")
    def capital_cost_regen_tank_constraint(b):
        regen_tank_vol_gal_dim = pyunits.convert(
            regen_tank_vol_gal * pyunits.gal**-1, to_units=pyunits.dimensionless
        )
        return (
            b.capital_cost_regen_tank
            == b.regen_tank_A_coeff * regen_tank_vol_gal_dim**b.regen_tank_b_coeff
        )

    @blk.Constraint(doc="Capital cost of backwashing/rinse tank")
    def capital_cost_backwash_tank_constraint(b):
        rinse_backwashing_vol_gal_dim = pyunits.convert(
            rinse_backwashing_vol_gal * pyunits.gal**-1,
            to_units=pyunits.dimensionless,
        )
        return (
            b.capital_cost_backwash_tank
            == b.backwash_tank_A_coeff
            * rinse_backwashing_vol_gal_dim**b.backwash_tank_b_coeff
        )

    @blk.Constraint(doc="Resin replacement operating cost")
    def resin_replacement_operating_cost_constraint(b):
        return b.resin_replacement_operating_cost == pyunits.convert(
            bed_vol_ft3
            * b.unit_model.number_columns_total
            * b.annual_resin_replacement_factor
            * 0.5
            * (b.anion_exchange_resin_cost + b.cation_exchange_resin_cost),
            to_units=b.basis_currency / pyunits.year,
        )

    @blk.Constraint(doc="Capital cost")
    def capital_cost_constraint(b):
        return b.capital_cost == pyunits.convert(
            (
                b.unit_model.number_ax_columns
                + b.unit_model.number_columns_redundant * 0.5
            )
            * (b.capital_cost_anion_exchange_resin + b.capital_cost_vessel),
            to_units=b.costing_package.base_currency,
        ) + pyunits.convert(
            (
                b.unit_model.number_cx_columns
                + b.unit_model.number_columns_redundant * 0.5
            )
            * (b.capital_cost_cation_exchange_resin + b.capital_cost_vessel),
            to_units=b.costing_package.base_currency,
        ) + pyunits.convert(
            b.capital_cost_regen_tank, to_units=b.costing_package.base_currency
        ) + pyunits.convert(
            b.capital_cost_backwash_tank, to_units=b.costing_package.base_currency
        )

    @blk.Constraint(doc="Fixed operating cost")
    def fixed_operating_cost_constraint(b):
        return b.fixed_operating_cost == pyunits.convert(
            b.resin_replacement_operating_cost,
            to_units=b.costing_package.base_currency / b.costing_package.base_period,
        )

    flow_mass_regen_solution = pyunits.convert(
        (
            (
                blk.unit_model.regeneration_dose
                * blk.unit_model.bed_volume
                * blk.unit_model.number_columns_total
            )
            / blk.unit_model.cycle_time
        )
        / blk.unit_model.regenerant_recycle,
        to_units=pyunits.kg / pyunits.s,
    )

    blk.cost_chemical_flow(
        material=blk.unit_model.regenerant, flow_expr=flow_mass_regen_solution
    )


@declare_process_block_class("IonExchange")
class UnitProcessData(WT3UnitProcessSIDOData):
    def build(self):
        super().build()


        if "target_component" in self.config.unit_params.keys():
            self.target_component = self.config.unit_params["target_component"]
            if self.target_component not in ["tds", "alkalinity", "alkalinity_as_caco3"]:
                raise ConfigurationError(f"'target_component' must be 'tds', 'alkalinity' or 'alkalinity_as_caco3' but {self.target_component} was provided.")
        else:
            self.target_component = "tds" # default is TDS
        
        target = self.target_component

        if "regenerant" in self.config.unit_params.keys():
            self.regenerant = self.config.unit_params["regenerant"]
        else:
            self.regenerant = "sodium_chloride" # default regenerant is NaCl

        # self.del_component(self.component_removal_equation)
        # self.del_component(self.component_remaining_equation)

        self.resin_capacity = Param(
            initialize=42,
            mutable=True, # typical operating capacity for resin is between 30-60 kg/m3; Wachowski, chap. 3
            units=pyunits.kg / pyunits.m**3,
            doc="Operational capacity of the resin in kg/m3 as CaCO3",
        )

        self.resin_density = Param(
            initialize=0.848,
            mutable=True,
            units=pyunits.kg / pyunits.m**3,
            doc="Resin density",
        )

        self.rinse_bed_volumes = Param(
            initialize=5,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Required number of bed volumes of rinse",
        )

        self.service_to_regen_flow_ratio = Param(
            initialize=3,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Ratio of service flow rate to regeneration flow rate",
        )

        self.backwashing_time = Param(
            initialize=10,
            mutable=True,
            units=pyunits.minutes,
            doc="Backwashing time",
        )

        self.backwashing_rate = Param(
            initialize=5,
            mutable=True,
            units=pyunits.m * pyunits.hr**-1,
            doc="Loading rate for backwashing cycle",
        )

        self.breakthrough_time = Param(
            initialize=24,
            mutable=True,
            units=pyunits.hr,
            doc="Breakthrough time (service cycle)",
        )

        self.regeneration_time = Param(
            initialize=0.5,
            mutable=True,
            units=pyunits.hr,
            doc="Time for one regeneration cycle",
        )

        self.regeneration_rate = Param(
            initialize=2,
            mutable=True,  # range is 0.25-0.4 gpm per cubic ft of resin - EPA-WBS
            units=pyunits.hr**-1,
            doc="Loading rate for regeneration cycle per cubic meter of resin",
        )

        self.regeneration_dose = Param(
            initialize=300,
            mutable=True,
            units=pyunits.kg / pyunits.m**3,
            doc="Mass of regenerant required per volume of resin",
        )

        self.regenerant_recycle = Param(
            initialize=1,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Number of times the regeneration solution can be recycled",
        )

        self.backwash_expansion_frac = Param(
            initialize=0.5,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Fraction of bed depth expanded during backwash",
        )

        self.ebct = Param(
            initialize=3,
            mutable=True,
            units=pyunits.min,
            doc="Empty bed conctact time",
        )

        self.free_board = Param(
            initialize=1,
            mutable=True,
            units=pyunits.m,
            doc="Free board height for distributor, underdrain",
        )

        self.redundant_column_freq = Param(
            initialize=4,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Frequency for redundant columns",
        )

        #################

        self.column_height = Var(
            initialize=1,
            bounds=(0, 4.26),  # EPA-WBS guidance
            units=pyunits.m,
            doc="Height of single column",
        )

        self.bed_depth = Var(
            initialize=1,
            bounds=(0.75, 2),  # EPA-WBS guidance
            units=pyunits.m,
            doc="Bed depth",
        )

        self.bed_diameter = Var(
            initialize=1,
            bounds=(0.75, 4.26),  # EPA-WBS guidance
            units=pyunits.m,
            doc="Bed diameter",
        )

        self.bed_volume_total = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Total bed volume required for each resin type",
        )

        self.bed_volume = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Bed volume per column",
        )

        self.loading_rate = Var(
            initialize=0.005,
            bounds=(0.0027, 0.011),
            units=pyunits.m / pyunits.s,
            doc="Loading rate (superficial bed velocity)",
        )

        self.number_ax_columns = Var(
            initialize=2,
            bounds=(1, None),
            units=pyunits.dimensionless,
            doc="Number of operational anion exchange columns",
        )

        self.number_cx_columns = Var(
            initialize=2,
            bounds=(1, None),
            units=pyunits.dimensionless,
            doc="Number of operational cation exchange columns",
        )

        self.number_columns_redundant = Var(
            initialize=1,
            bounds=(1, None),
            units=pyunits.dimensionless,
            doc="Number of redundant columns for each type.",
        )

        self.handle_unit_params()

        @self.Expression(doc="Total number of operational columns")
        def number_columns_op(b):
            return b.number_cx_columns + b.number_ax_columns

        @self.Expression(doc="Total number of columns")
        def number_columns_total(b):
            return b.number_columns_op + b.number_columns_redundant

        @self.Expression(doc="Mass of target to be removed per cycle")
        def mass_removed_per_cycle(b):
            return pyunits.convert(
                b.properties_out.flow_mass_comp[target] * b.breakthrough_time,
                to_units=pyunits.kg,
            )

        @self.Expression(doc="Cross-sectional area of one column")
        def bed_area(b):
            return pyunits.convert(
                Constants.pi * (b.bed_diameter / 2) ** 2, to_units=pyunits.m**2
            )

        @self.Expression(doc="Regeneration flow rate")
        def regen_flow_rate(b):
            return pyunits.convert(
                b.bed_volume_total * b.regeneration_rate,
                to_units=pyunits.m**3 / pyunits.s,
            )

        @self.Expression(doc="Regeneration volume required per cycle")
        def regeneration_volume(b):
            return pyunits.convert(
                b.regen_flow_rate * b.regeneration_time, to_units=pyunits.m**3
            )

        @self.Expression(doc="Time required for rinsing")
        def rinse_time(b):
            return pyunits.convert(b.ebct * b.rinse_bed_volumes, to_units=pyunits.min)

        @self.Expression(doc="Rinse flow rate")
        def rinse_rate(b):
            return pyunits.convert(
                b.properties_in.flow_vol / b.service_to_regen_flow_ratio,
                to_units=pyunits.m**3 / pyunits.s,
            )

        @self.Expression(doc="Rinse volume required per cycle")
        def rinse_volume(b):
            return pyunits.convert(b.rinse_rate * b.rinse_time, to_units=pyunits.m**3)

        @self.Expression(doc="Backwashing volume required per cycle")
        def backwashing_volume(b):
            return (
                pyunits.convert(
                    b.backwashing_rate * b.bed_area * b.backwashing_time,
                    to_units=pyunits.m**3,
                )
                * b.number_columns_op
            )

        @self.Expression(doc="Total cycle time")
        def cycle_time(b):
            return (
                pyunits.convert(b.breakthrough_time, to_units=pyunits.s)
                + pyunits.convert(b.regeneration_time, to_units=pyunits.s)
                + pyunits.convert(b.rinse_time, to_units=pyunits.s)
                + pyunits.convert(b.backwashing_time, to_units=pyunits.s)
            )

        @self.Expression(doc="Column volume")
        def column_volume(b):
            return pyunits.convert(
                b.bed_area * b.column_height, to_units=pyunits.m**3
            )

        @self.Constraint(doc="Total bed volume required")
        def bed_volume_total_constraint(b):
            return b.bed_volume_total >= pyunits.convert(
                b.mass_removed_per_cycle / b.resin_capacity, to_units=pyunits.m**3
            )

        @self.Constraint(doc="Bed volume per operational column")
        def bed_volume_constraint(b):
            return b.bed_volume == b.bed_area * b.bed_depth

        @self.Constraint(doc="Number of operational columns")
        def number_operational_columns_constraint(b):
            return b.number_ax_columns * b.bed_volume == b.bed_volume_total

        @self.Constraint(
            doc="Number of anion columns is equal to the number of cation columns"
        )
        def number_cation_anion_columns_constraint(b):
            return b.number_ax_columns == b.number_cx_columns

        @self.Constraint(doc="Number of redundant columns")
        def number_redundant_columns_constraint(b):
            return (
                b.number_columns_redundant
                >= b.number_columns_op / b.redundant_column_freq
            )

        @self.Constraint(doc="Loading rate")
        def loading_rate_constraint(b):
            return b.loading_rate == pyunits.convert(
                b.bed_depth / b.ebct, to_units=pyunits.m / pyunits.s
            )

        @self.Constraint(doc="Height of column")
        def column_height_constraint(b):
            return (
                b.column_height
                == b.bed_depth * (1 + b.backwash_expansion_frac) + b.free_board
            )

        @self.Constraint(doc="Flow per column")
        def flow_per_column_constraint(b):
            return b.properties_in.flow_vol <= pyunits.convert(
                b.loading_rate * b.bed_area * b.number_ax_columns,
                to_units=pyunits.m**3 / pyunits.s,
            )

        flow_mass_regen_solution = pyunits.convert(
            (
                (self.regeneration_dose * self.bed_volume * self.number_columns_total)
                / self.cycle_time
            )
            / self.regenerant_recycle,
            to_units=pyunits.kg / pyunits.s,
        )

        # @self.Constraint(
        #     self.config.property_package.solute_set,
        #     doc="Solute removal equation for ion exchange",
        # )
        # def solute_removal_equation(b, j):
        #     if j == target:
        #         return (
        #             b.removal_fraction[j] * b.properties_in.flow_mass_comp[j]
        #             + flow_mass_regen_solution
        #             == b.properties_waste.flow_mass_comp[j]
        #         )
        #     else:
        #         return (
        #             b.removal_fraction[j] * b.properties_in.flow_mass_comp[j]
        #             == b.properties_waste.flow_mass_comp[j]
        #         )

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
            if k == "flow_mass_comp":
                # elif isinstance(v, dict):
                state_args_out[k] == dict()
                for j, u in v.items():
                    if j == "H2O":
                        state_args_out[k][j] = u
                    else:
                        state_args_out[k][j] = (1 - self.removal_fraction[j].value) * u

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
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

        # self.bed_volume_total.fix()
        ncol = math.ceil(value(self.number_columns_op))
        ncol_red = ncol / value(self.redundant_column_freq)
        self.number_ax_columns.set_value(ncol)
        self.number_columns_redundant.set_value(ncol_red)
        # self.bed_depth.fix()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)

        init_log.info("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        
        # Release Inlet state
        self.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        # if not check_optimal_termination(res):
        #     raise InitializationError(f"Unit model {self.name} failed to initialize.")

        self.initialized = True

    @property
    def default_costing_method(self):
        return cost_ion_exchange
