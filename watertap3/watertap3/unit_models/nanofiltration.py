from pyomo.environ import Var, Param, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from watertap3.core.wt3_unit_sido import WT3UnitProcessSIDOData

## REFERENCE
## CAPITAL:
# Texas Water Board User's Manual for Integrated Treatment Train Toolbox - Potable Reuse (IT3PR) Version 2.0.
## ELECTRICITY:
# Plappally, A. K., & Lienhard V, J. H. (2012). doi:10.1016/j.rser.2012.05.022

module_name = "nanofiltration"


def cost_nanofiltration(blk):
    blk.basis_year = 2014
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")

    blk.capital_cost_base = Var(
        initialize=2500000,
        bounds=(0, None),
        units=blk.basis_currency * pyunits.day * pyunits.Mgallons**-1,
        doc="Nanofiltration capital cost equation slope",
    )

    blk.handle_costing_unit_params()
    blk.fix_all_vars()
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")

    @blk.Constraint(doc="Capital cost constraint")
    def capital_cost_constraint(b):
        flow_mgd = pyunits.convert(
            b.unit_model.properties_in.flow_vol, to_units=pyunits.Mgallons / pyunits.day
        )
        return b.capital_cost == pyunits.convert(
            b.capital_cost_base * flow_mgd, to_units=b.costing_package.base_currency
        )

    @blk.Expression(doc="Pumping power required")
    def power_required(b):
        return (
            pyunits.convert(
                b.unit_model.properties_in.flow_vol * b.unit_model.operating_pressure,
                to_units=pyunits.kilowatt,
            )
            / b.unit_model.pump_efficiency
        )

    blk.costing_package.cost_flow(blk.power_required, "electricity")


@declare_process_block_class("Nanofiltration")
class UnitProcessData(WT3UnitProcessSIDOData):
    def build(self):
        super().build()

        self.operating_pressure = Param(
            initialize=10,
            mutable=True,
            units=pyunits.bar,
            doc="Nanofiltration operating pressure",
        )

        self.pump_efficiency = Param(
            initialize=0.8,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Nanofiltration pump efficiency",
        )

        self.handle_unit_params()

    @property
    def default_costing_method(self):
        return cost_nanofiltration
