from pyomo.environ import Expression, units as pyunits
from watertap3.utils import epa_cost_curve, financials
from watertap3.core.wt3_unit_sido import WT3UnitProcess

## REFERENCE: ADD REFERENCE HERE
# CAPITAL AND ELECTRICITY:
# EPA Drinking Water Treatment Technology Unit Cost Models
# https://www.epa.gov/sdwa/drinking-water-treatment-technology-unit-cost-models
# EPA model run several times under different conditions
# Regression for capital and electricity costs based on regression from those model outputs.

module_name = 'packed_tower_aeration'

class UnitProcess(WT3UnitProcess):

    def fixed_cap(self):
        '''

        :return: Fixed capital for packed tower aeration [$MM]
        '''
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.hr)
        if 'radon_rem' in self.unit_params.keys():
            self.radon_rem = self.unit_params['radon_rem']
        else:
            self.radon_rem = 0.99
        self.cost_coeffs, self.elect_coeffs, self.mats_name, self.mats_cost, _ = \
            epa_cost_curve(module_name, radon_rem=self.radon_rem)
        for k, v in self.mats_cost.items():
            self.mats_cost[k] = v * (pyunits.kg/pyunits.m**3)
        self.chem_dict = self.mats_cost
        pta_cap = self.cost_coeffs[0] * self.flow_in ** self.cost_coeffs[1]
        return pta_cap * 1E-6

    def elect(self):
        '''
        Electricity intensity for packed tower aeration [kWh/m3]
        :return:
        '''
        electricity = self.elect_coeffs[0] * self.flow_in ** self.elect_coeffs[1]
        return electricity

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2017
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.elect(),
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year)