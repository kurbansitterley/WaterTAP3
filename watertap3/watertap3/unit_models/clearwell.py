from copy import deepcopy

from pyomo.environ import (
    Param,
    Var,
    check_optimal_termination,
    value,
    units as pyunits,
)

import idaes.logger as idaeslog
from idaes.core import declare_process_block_class
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.scaling import set_scaling_factor, get_scaling_factor

from watertap.costing.util import make_capital_cost_var

from watertap3.core.wt3_unit_pt import WT3UnitProcessPTData


## REFERENCE
## CAPITAL:
# Kawamura and McGivney
## ELECTRICITY:
# citation here

__author__ = "Kurban Sitterley"
module_name = "clearwell"


def cost_clearwell(blk):

    blk.basis_year = 2007
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")

    blk.capital_cost_base = Var(
        initialize=1016775,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Clearwell capital cost base",
    )

    blk.capital_cost_exp = Var(
        initialize=0.6836,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Clearwell capital cost exponent",
    )

    blk.fix_all_vars()

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)

    @blk.Constraint(doc="Capital cost for clearwell")
    def capital_cost_constraint(b):
        storage_vol_dim = pyunits.convert(
            b.unit_model.storage_volume * pyunits.Mgallons**-1,
            to_units=pyunits.dimensionless,
        )

        return b.capital_cost == pyunits.convert(
            b.capital_cost_base * storage_vol_dim**b.capital_cost_exp,
            to_units=blk.costing_package.basis_currency,
        )


@declare_process_block_class("Clearwell")
class UnitProcessData(WT3UnitProcessPTData):
    def build(self):
        super().build()

        self.storage_time = Param(
            initialize=3,
            mutable=True,
            units=pyunits.hour,
            doc="Storage time for clearwell",
        )

        self.storage_volume = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.Mgallons,
            doc="Storage volume for clearwell",
        )

        @self.Constraint(doc="Storage volume calculation")
        def storage_volume_constraint(b):
            return b.storage_volume == pyunits.convert(
                b.properties_in.flow_vol * b.storage_time, to_units=pyunits.Mgallons
            )

        self.handle_unit_params()

    @property
    def default_costing_method(self):
        return cost_clearwell
