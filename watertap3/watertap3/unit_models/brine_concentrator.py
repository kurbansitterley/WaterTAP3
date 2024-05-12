from pyomo.environ import Var, Constraint, Param, Expression, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from watertap3.core.wt3_unit_sido import WT3UnitProcessSIDOData
from watertap3.core.util.pumping_energy import pumping_energy

## REFERENCE
## CAPITAL AND ELECTRICITY:
# Survey of High-Recovery and Zero Liquid Discharge Technologies for Water Utilities (2008).
# WateReuse Foundation
# https://www.waterboards.ca.gov/water_issues/programs/grants_loans/water_recycling/research/02_006a_01.pdf
# Regressions for capital and electricity developed from
# data in Table 5.1, Table A2.3
# Capital = f(TDS, recovery, flow)
# Electricity = f(TDS, recovery, flow)

module_name = "brine_concentrator"


def cost_brine_concentrator(blk):
    blk.basis_year = 2007
    blk.basis_currency = getattr(pyunits, f"MUSD_{blk.basis_year}")

    tds_in = pyunits.convert(
        blk.unit_model.properties_in.conc_mass_comp["tds"],
        to_units=pyunits.mg / pyunits.liter,
    )
    flow_in = pyunits.convert(
        blk.unit_model.properties_in.flow_vol, to_units=pyunits.m**3 / pyunits.hr
    )

    # brine concetrator capital = 15.1 + (tds_in * 3.02E-4) - (water_recovery * 18.8) + (flow_in * 8.08E-2)

    blk.capital_cost_intercept = Var(
        initialize=15.1,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Intercept for capital cost equation",
    )
    blk.capital_cost_tds_coeff = Var(
        initialize=3.02e-4,
        bounds=(0, None),
        units=blk.basis_currency * pyunits.liter * pyunits.mg**-1,
        doc="TDS coefficient for capital cost equation",
    )
    blk.capital_cost_water_recovery_coeff = Var(
        initialize=-18.8,
        bounds=(None, None),
        units=blk.basis_currency,
        doc="Water recovery coefficient for capital cost equation",
    )
    blk.capital_cost_flow_coeff = Var(
        initialize=8.08e-2,
        bounds=(0, None),
        units=blk.basis_currency * pyunits.hr * pyunits.m**-3,
        doc="Volumetric flow coefficient for capital cost equation",
    )

    # brine concentrator energy intensity = 9.73 + tds_in * 1.1E-4 + water_recovery * 10.4 + flow_in * 3.83E-5

    blk.energy_intensity_intercept = Var(
        initialize=9.73,
        bounds=(0, None),
        units=pyunits.kWh * pyunits.m**-3,
        doc="Intercept for energy intensity equation",
    )
    blk.energy_intensity_tds_coeff = Var(
        initialize=1.1e-4,
        bounds=(0, None),
        units=pyunits.kWh * pyunits.m**-3 * pyunits.liter * pyunits.mg**-1,
        doc="TDS coefficient for energy intensity equation",
    )
    blk.energy_intensity_water_recovery_coeff = Var(
        initialize=10.4,
        bounds=(0, None),
        units=pyunits.kWh * pyunits.m**-3,
        doc="Water recovery coefficient for energy intensity equation",
    )
    blk.energy_intensity_flow_coeff = Var(
        initialize=3.83e-5,
        bounds=(0, None),
        units=pyunits.kWh * pyunits.hr * pyunits.m**-6,
        doc="Volumetric flow coefficient for energy intensity equation",
    )

    blk.handle_costing_unit_params()
    blk.fix_all_vars()
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)

    blk.energy_intensity = Var(
        initialize=0.1,
        bounds=(0, None),
        units=pyunits.kWh / pyunits.m**3,
        doc="Chlorination energy intensity",
    )

    @blk.Constraint(doc="Capital cost equation for brine concentrator")
    def capital_cost_constraint(b):
        return blk.capital_cost == pyunits.convert(
            b.capital_cost_intercept
            + b.capital_cost_tds_coeff * tds_in
            + b.capital_cost_water_recovery_coeff * b.unit_model.water_recovery
            + b.capital_cost_flow_coeff * flow_in,
            to_units=blk.costing_package.base_currency,
        )

    @blk.Constraint(doc="Energy intensity for brine concentrator")
    def energy_intensity_constraint(b):
        return blk.energy_intensity == pyunits.convert(
            b.energy_intensity_intercept
            + b.energy_intensity_tds_coeff * tds_in
            + b.energy_intensity_water_recovery_coeff * b.unit_model.water_recovery
            + b.energy_intensity_flow_coeff * flow_in,
            to_units=pyunits.kWh / pyunits.m**3,
        )

    @blk.Expression(doc="Power required")
    def power_required(b):
        return pyunits.convert(
            b.energy_intensity * b.unit_model.properties_in.flow_vol,
            to_units=pyunits.kilowatt,
        )

    blk.costing_package.cost_flow(blk.power_required, "electricity")


@declare_process_block_class("BrineConcentrator")
class UnitProcessData(WT3UnitProcessSIDOData):
    def build(self):
        super().build()

    @property
    def default_costing_method(self):
        return cost_brine_concentrator
