from pyomo.environ import Expression
from watertap3.utils import financials
from watertap3.core.wt3_unit_pt import WT3UnitProcessPTData
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var, make_fixed_operating_cost_var, register_costing_parameter_block

module_name = 'passthrough'

def cost_passthrough(blk):

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)
    blk.capital_cost.fix(0)


@declare_process_block_class("UnitProcess")
class UnitProcessData(WT3UnitProcessPTData):

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=0,
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=0,
                doc='Electricity intensity [kWh/m3]')
        self.deltaP_outlet.unfix()
        # self.deltaP_waste.unfix()
        self.pressure_out.fix(1)
        # self.pressure_waste.fix(1)
        financials.get_complete_costing(self.costing)

    @property
    def default_costing_method(self):
        return cost_passthrough