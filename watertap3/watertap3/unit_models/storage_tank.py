from pyomo.environ import Var, Param, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from watertap3.core.wt3_unit_pt import WT3UnitProcessPTData


## REFERENCE
## CAPITAL:
# "Cone roof tank" costs from:
# DOE/NETL-2002/1169 - Process Equipment Cost Estimation Final Report
# Loh, H. P., Lyons, Jennifer, and White, Charles W. Process Equipment Cost Estimation, Final Report.
# United States: N. p., 2002. Web. doi:10.2172/797810.
# Regression of cost vs. capacity
# Capacity calculated based on storage time (user input)

module_name = "storage_tank"


def cost_storage_tank(blk):
    blk.basis_year = 1998
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")

    blk.capital_cost_base = Var(
        initialize=0.00344,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Storage tank capital cost basis",
    )

    blk.capital_cost_exp = Var(
        initialize=0.72093,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Storage tank capital cost exponent",
    )

    blk.handle_costing_unit_params()

    blk.fix_all_vars()

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)

    @blk.Constraint(doc="Capital cost equation")
    def capital_cost_constraint(b):
        tank_vol_dim = pyunits.convert(
            blk.unit_model.storage_volume * pyunits.m**-3,
            to_units=pyunits.dimensionless,
        )
        return b.capital_cost == pyunits.convert(
            b.capital_cost_base * tank_vol_dim**b.capital_cost_exp,
            to_units=b.costing_package.base_currency,
        )


@declare_process_block_class("StorageTank")
class UnitProcessData(WT3UnitProcessPTData):
    def build(self):
        super().build()

        self.storage_time = Param(
            initialize=6,
            mutable=True,
            units=pyunits.hour,
            doc="Storage duration",
        )

        self.surge_capacity = Param(
            initialize=0.2,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Storage tank surge capacity",
        )

        self.handle_unit_params()

        self.storage_volume = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Storage tank volume",
        )

        @self.Constraint(doc="Storage volume equation")
        def storage_volume_constraint(b):
            return b.storage_volume == pyunits.convert(
                b.properties_in.flow_vol * b.storage_time * (1 + b.surge_capacity),
                to_units=pyunits.m**3,
            )

    @property
    def default_costing_method(self):
        return cost_storage_tank
