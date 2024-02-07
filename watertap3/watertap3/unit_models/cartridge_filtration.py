from pyomo.environ import Expression, units as pyunits
from watertap3.utils import financials
from watertap3.core.wt3_unit import WT3UnitProcess

## REFERENCE:
# from IT3PR, section 3.5.6 figure 3.3

module_name = 'cartridge_filtration'

class UnitProcess(WT3UnitProcess):

    def fixed_cap(self):
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.Mgallons / pyunits.day)
        
        # self.base_fixed_cap_cost = 0.72557
        # self.cap_scaling_exp = 0.5862
        self.base_fixed_cap_cost = 0.2161
        self.cap_scaling_exp = 0.7639
        cart_filt_cap = self.base_fixed_cap_cost * self.flow_in ** self.cap_scaling_exp
        return cart_filt_cap

    def elect(self):
        electricity = 2E-4
        return electricity

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2014
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.elect(),
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year)