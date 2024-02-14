from watertap3.core.wt3_unit_pt import WT3UnitProcessPTData
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var

module_name = 'passthrough'

def cost_passthrough(blk):

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)
    blk.capital_cost.fix(0)


@declare_process_block_class("Passthrough")
class UnitProcessData(WT3UnitProcessPTData):

    def build(self):
        super().build()

    @property
    def default_costing_method(self):
        return cost_passthrough