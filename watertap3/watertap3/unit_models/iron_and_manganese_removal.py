from pyomo.environ import Var, Param, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from idaes.core.util.scaling import set_scaling_factor, get_scaling_factor
from watertap3.core.wt3_unit_siso import WT3UnitProcessSISOData

## REFERENCE
## CAPITAL:
# Based on costs for FILTER MEDIA DUAL MEDIA - FIGURE 5.5.27
# Cost Estimating Manual for Water Treatment Facilities (McGivney/Kawamura) (2008)
# DOI:10.1002/9780470260036
# https://www.lenntech.com/schema-of-an-iron-removal-system.htm
## ELECTRICITY:
# FOR BLOWER
# Loh, H. P., Lyons, Jennifer, and White, Charles W. Process Equipment Cost Estimation, Final Report.
# United States: N. p., 2002. Web. doi:10.2172/797810.

module_name = "iron_and_manganese_removal"


def cost_iron_and_manganese_removal(blk):

    blk.basis_year = 2014
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")

    blk.flow_scaling_base = Var(
        initialize=4732,
        bounds=(0, None),
        units=pyunits.m**3 / pyunits.hr,
        doc="Flow scaling base",
    )
    blk.capital_cost_exp = Var(
        initialize=0.7,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Fe/Mn removal capital cost exponent",
    )
    blk.filter_capital_cost_intercept = Var(
        initialize=21377,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Dual media filtration capital cost equation - intercept",
    )
    blk.filter_capital_cost_slope = Var(
        initialize=38.319,
        bounds=(0, None),
        units=blk.basis_currency / pyunits.ft**2,
        doc="Dual media filtration capital cost equation - slope",
    )
    blk.backwash_capital_cost_intercept = Var(
        initialize=92947,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Backwashing system capital cost equation - intercept",
    )
    blk.backwash_capital_cost_slope = Var(
        initialize=292.44,
        bounds=(0, None),
        units=blk.basis_currency / pyunits.ft**2,
        doc="Backwashing system capital cost equation - slope",
    )
    blk.energy_intensity_blower = Var(
        initialize=110.2,  # converted from 147.8 hp / (m3/hr)
        bounds=(0, None),
        units=pyunits.kWh / pyunits.m**3,
        doc="Energy intensity of blower",
    )
    blk.blower_capital_cost = Var(
        initialize=100000,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Blower capital cost",
    )

    blk.handle_costing_unit_params()

    blk.fix_all_vars()

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TPEC")
    filter_area_ft2 = pyunits.convert(
        blk.unit_model.filter_surface_area, to_units=pyunits.ft**2
    )

    blk.capital_cost_base = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Fe/Mn capital cost base",
    )
    blk.filter_capital_cost = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Filtration capital cost",
    )
    blk.backwash_capital_cost = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Backwash system capital cost",
    )

    @blk.Constraint(doc="Filter capital cost equation")
    def filter_capital_cost_constraint(b):
        return (
            b.filter_capital_cost
            == (
                b.filter_capital_cost_intercept
                + b.filter_capital_cost_slope * filter_area_ft2
            )
            * b.unit_model.number_units
        )

    @blk.Constraint(doc="Backwashing system capital cost equation")
    def backwash_capital_cost_constraint(b):
        return (
            b.backwash_capital_cost
            == b.backwash_capital_cost_intercept
            + b.backwash_capital_cost_slope * filter_area_ft2
        )

    @blk.Constraint(doc="Base capital cost")
    def capital_cost_base_constraint(b):
        return b.capital_cost_base == b.costing_package.TPEC * (
            b.blower_capital_cost + b.filter_capital_cost + b.backwash_capital_cost
        )

    @blk.Constraint(doc="Unit total capital cost")
    def capital_cost_constraint(b):
        flow_m3_hr = pyunits.convert(
            b.unit_model.properties_in.flow_vol, to_units=pyunits.m**3 / pyunits.hr
        )
        return b.capital_cost == pyunits.convert(
            b.capital_cost_base
            * (flow_m3_hr / b.flow_scaling_base) ** b.capital_cost_exp,
            to_units=b.costing_package.base_currency,
        )

    @blk.Expression(doc="Blower power required")
    def power_required(b):
        return pyunits.convert(
            b.energy_intensity_blower * b.unit_model.air_flow_rate, to_units=pyunits.kW
        )

    blk.costing_package.cost_flow(blk.power_required, "electricity")


@declare_process_block_class("IronAndManganeseRemoval")
class UnitProcessData(WT3UnitProcessSISOData):
    def build(self):
        super().build()

        self.number_units = Param(
            initialize=6,
            units=pyunits.dimensionless,
            doc="Number of units",
        )
        self.filter_surface_area = Param(
            initialize=580,
            units=pyunits.m**2,
            doc="Filter surface area of unit",
        )
        self.air_water_ratio = Param(
            initialize=0.001,
            units=pyunits.dimensionless,
            doc="Air-to-water ratio for blower",
        )

        @self.Expression(doc="Air flow rate equation")
        def air_flow_rate(b):
            return b.air_water_ratio * b.properties_in.flow_vol

        self.handle_unit_params()

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if get_scaling_factor(self.air_flow_rate) is None:
            set_scaling_factor(self.air_flow_rate, 1e4)

    @property
    def default_costing_method(self):
        return cost_iron_and_manganese_removal
