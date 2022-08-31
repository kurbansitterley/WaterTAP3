from pyomo.environ import Expression, value, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# "Precious Metal Heap Leach Design and Practice" Daniel W. Kappes (2002)
# Mineral Processing Plant Design, Practice, and Control, Volume 1, Published 2002 by SME (1606)
# http://ore-max.com/pdfs/resources/precious_metal_heap_leach_design_and_practice.pdf

module_name = 'solution_distribution_and_recovery_plant'

class UnitProcess(WT3UnitProcess):

    def fixed_cap(self):
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.hr)
        try:
            self.mining_capacity = self.unit_params['mining_capacity'] * (pyunits.tonnes/pyunits.day)
            self.ore_heap_soln = self.unit_params['ore_heap_soln'] * (pyunits.gallons/pyunits.tonnes)
            self.make_up_water = 85 / 500 * self.ore_heap_soln * (pyunits.gallons/pyunits.tonnes)
            self.recycle_water = self.ore_heap_soln - self.make_up_water
            self.make_up_water = pyunits.convert(self.make_up_water * self.mining_capacity,
                to_units=(pyunits.m**3/pyunits.hour))
            self.recycle_water = pyunits.convert(self.recycle_water * self.mining_capacity,
                to_units=(pyunits.m**3/pyunits.hour))

        except:
            self.mining_capacity = 922 * (pyunits.tonnes/pyunits.day)
            self.ore_heap_soln = 500 * (pyunits.gallons/pyunits.tonnes)
            self.make_up_water = 85 * (pyunits.gallons/pyunits.tonnes)
            self.recycle_water = self.ore_heap_soln - self.make_up_water
            self.make_up_water = pyunits.convert(self.make_up_water * self.mining_capacity,
                to_units=(pyunits.m**3/pyunits.hour))
            self.recycle_water = pyunits.convert(self.recycle_water * self.mining_capacity,
                to_units=(pyunits.m**3/pyunits.hour))

        self.dist_recov = 0.00347 * self.mining_capacity ** 0.71917

        self.perc = 0.3102 * self.mining_capacity ** 0.1119  # regression made by KAS in excel - mining capacity vs percent of subtotal
        if value(self.perc) > 1:
            self.perc = 1

        self.dist_recov_basis = self.dist_recov * (1 + self.perc)

        self.dist_recov_op = 7.71759 * self.mining_capacity ** 0.91475
        self.dist_recov_exp = 0.719
        self.dist_recov_other = self.dist_recov_op * 1E-6 * 365
        self.costing.other_var_cost = self.dist_recov_other

        self.flow_factor = self.flow_in / self.recycle_water
        
        dist_recov_cap = self.flow_factor * self.dist_recov_basis ** self.dist_recov_exp
        return dist_recov_cap

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