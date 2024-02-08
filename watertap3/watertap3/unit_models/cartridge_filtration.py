from pyomo.environ import Var, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from watertap3.core.wt3_unit_siso import WT3UnitProcessSISOData
## REFERENCE:
# from IT3PR, section 3.5.6 figure 3.3

module_name = "cartridge_filtration"


@declare_process_block_class("UnitProcess")
class UnitProcessData(WT3UnitProcessSISOData):
    def build(self):
        super().build()

    @property
    def default_costing_method(self):
        return cost_cartridge_filtration


def cost_cartridge_filtration(blk):

    unit_params = blk.unit_model.unit_params

    blk.basis_year = 2014
    blk.basis_currency = getattr(pyunits, f"MUSD_{blk.basis_year}")

    blk.capital_cost_base = Var(
        initialize=0.2161,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Cartridge filtration capital cost basis",
    )
    blk.capital_cost_exp = Var(
        initialize=0.7639,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Cartridge filtration capital cost exponent",
    )

    blk.energy_intensity = Var(
        initialize=2e-4,
        bounds=(0, None),
        units=pyunits.kWh / pyunits.m**3,
        doc="Cartridge filtration energy intensity",
    )

    blk.fix_all_vars()

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)

    flow_mgd = pyunits.convert(
        blk.unit_model.properties_in.flow_vol, to_units=pyunits.Mgallons / pyunits.day
    )

    @blk.Constraint(doc="Unit total capital cost")
    def captital_cost_constraint(b):
        flow_dim = pyunits.convert(
            flow_mgd * pyunits.day / pyunits.Mgallons, to_units=pyunits.dimensionless
        )
        return b.capital_cost == pyunits.convert(
            b.capital_cost_base * flow_dim**b.capital_cost_exp,
            to_units=b.costing_package.base_currency,
        )

    @blk.Expression(doc="Power required")
    def power_required(b):
        return pyunits.convert(
            b.energy_intensity * b.unit_model.properties_in.flow_vol,
            to_units=pyunits.kilowatt,
        )

    blk.costing_package.cost_flow(blk.power_required, "electricity")
