from pyomo.environ import Var, Constraint, Param, Expression, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from watertap3.core.wt3_unit_pt import WT3UnitProcessPTData
from watertap3.core.util.pumping_energy import pumping_energy

## REFERENCE
## CAPITAL:
# citation here
## ELECTRICITY:
# citation here

module_name = 'unit'

def cost_unit_model(blk):
    blk.basis_year = 2018
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")
    blk.capital_cost_base = Var(
        initialize=123,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Unit capital cost basis",
    )
    blk.capital_cost_exp = Var(
        initialize=0.7,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Unit capital cost exponent",
    )
    ## VARS GO HERE
    blk.handle_costing_unit_params()
    blk.fix_all_vars()
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)
    ## CONSTRAINTS GO HERE

@declare_process_block_class("UnitProcess")
class UnitProcessData(WT3UnitProcessPTData):

    def build(self):
        super().build()
    
    @property
    def default_costing_method(self):
        return cost_unit_model