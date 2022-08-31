from pyomo.environ import Expression
from watertap3.utils import financials
from watertap3.wt_units.wt_unit_pt import WT3UnitProcessPT

module_name = 'passthrough'

class UnitProcess(WT3UnitProcessPT):

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