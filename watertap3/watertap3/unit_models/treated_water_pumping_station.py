from pyomo.environ import Var, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from watertap3.core.wt3_unit_pt import WT3UnitProcessPTData

## REFERENCE: McGivney & Kawamura (2008) - Figure 5.5.35b

module_name = "treated_water_pumping_station"


def cost_treated_water_pumping_station(blk):
    blk.basis_year = 2007
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")

    blk.capital_cost_base = Var(
        initialize=45314,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Municipal drinking capital cost basis",
    )

    blk.capital_cost_exp = Var(
        initialize=0.84,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Municipal drinking capital cost exponent",
    )

    blk.add_pumping_energy(lift_height=300)
    blk.handle_costing_unit_params()

    blk.fix_all_vars()

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TPEC")
    flow_in_mgd = pyunits.convert(
        blk.unit_model.properties_in.flow_vol, to_units=pyunits.Mgallons / pyunits.day
    )

    @blk.Constraint(doc="Capital cost equation")
    def capital_cost_constraint(b):
        flow_in_dim = pyunits.convert(
            flow_in_mgd * pyunits.day * pyunits.Mgallons**-1,
            to_units=pyunits.dimensionless,
        )
        return b.capital_cost == b.costing_package.TPEC * pyunits.convert(
            b.capital_cost_base * flow_in_dim**b.capital_cost_exp,
            to_units=blk.costing_package.base_currency,
        )


@declare_process_block_class("TreatedWaterPumpingStation")
class UnitProcessData(WT3UnitProcessPTData):
    def build(self):
        super().build()

    @property
    def default_costing_method(self):
        return cost_treated_water_pumping_station
