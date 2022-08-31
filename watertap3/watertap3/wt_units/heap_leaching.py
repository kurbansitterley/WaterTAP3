from pyomo.environ import Expression, value, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# "Precious Metal Heap Leach Design and Practice" Daniel W. Kappes (2002)
# Mineral Processing Plant Design, Practice, and Control, Volume 1, Published 2002 by SME (1606)
# http://ore-max.com/pdfs/resources/precious_metal_heap_leach_design_and_practice.pdf


module_name = 'heap_leaching'

class UnitProcess(WT3UnitProcess):

    def fixed_cap(self):
        time = self.flowsheet().config.time.first()
        
        self.flow_in = pyunits.convert(self.flow_vol_in[time], to_units=pyunits.m**3/pyunits.hr)
        try:
            self.mining_capacity = self.unit_params['mining_capacity'] * (pyunits.tonnes/pyunits.day)
            self.ore_heap_soln = self.unit_params['ore_heap_soln'] * (pyunits.gallons/pyunits.tonnes)
            self.make_up_water = 85 / 500 * self.ore_heap_soln * (pyunits.gallons/pyunits.tonnes)
            self.make_up_water = pyunits.convert(self.make_up_water * self.mining_capacity,
                to_units=(pyunits.m**3/pyunits.hour))
        except:
            self.mining_capacity = 922 * (pyunits.tonnes/pyunits.day)
            self.ore_heap_soln = 500 * (pyunits.gallons/pyunits.tonnes)
            self.make_up_water = 85 * (pyunits.gallons/pyunits.tonnes)
            self.make_up_water = pyunits.convert(self.make_up_water * self.mining_capacity,
                to_units=(pyunits.m**3/pyunits.hour))
        self.mine_equip = 0.00124 * self.mining_capacity ** 0.93454
        self.mine_develop = 0.01908 * self.mining_capacity ** 0.43068
        self.crushing = 0.0058 * self.mining_capacity ** 0.6651
        self.leach = 0.0005 * self.mining_capacity ** 0.94819
        self.stacking = 0.00197 * self.mining_capacity ** 0.77839
        self.dist_recov = 0.00347 * self.mining_capacity ** 0.71917
        self.subtotal = (self.mine_equip + self.mine_develop + self.crushing + \
            self.leach + self.stacking + self.dist_recov)
        self.perc = 0.3102 * self.mining_capacity ** 0.1119  # regression made by KAS in excel - mining capacity vs percent of subtotal
        if value(self.perc) > 1:
            self.perc = 1

        self.mining_to_heap_basis = (self.mine_equip + self.mine_develop + self.crushing + \
            self.leach) * (1 + self.perc)
        self.mining_to_heap_exp = (self.mine_equip * 0.935 + self.mine_develop * 0.431 + \
            self.crushing * 0.665 + self.leach * 0.948) / \
                (self.mine_equip + self.mine_develop + self.crushing + self.leach)

        self.mine_equip_op = 22.54816 * self.mining_capacity ** 0.74807
        self.crushing_op = 4.466 * self.mining_capacity ** 0.8794
        self.leach_op = 6.34727 * self.mining_capacity ** 0.68261

        # self.mining_to_heap_other = (self.mine_equip_op + self.crushing_op + self.leach_op) / self.make_up_water
        self.mining_to_heap_other = (self.mine_equip_op + self.crushing_op + self.leach_op) * 1E-6 * 365 # cost curves return $/day
        self.costing.other_var_cost = self.mining_to_heap_other
        self.flow_factor = self.flow_in / 73
        heap_cap = self.flow_factor * self.mining_to_heap_basis ** self.mining_to_heap_exp
        return heap_cap

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2008
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment') 
        self.electricity = Expression(expr=0,
                doc='Electricity intensity [kWh/m3]')  
        financials.get_complete_costing(self.costing, basis_year=basis_year)


