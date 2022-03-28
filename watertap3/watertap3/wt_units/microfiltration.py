from pyomo.environ import Expression, Var, Constraint, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE: TWB

module_name = 'microfiltration'
tpec_or_tic = 'TPEC'

class UnitProcess(WT3UnitProcess):

    def fixed_cap(self):
        time = self.flowsheet().config.time.first()
        if self.cost_method == 'twb':
            self.chem_dict = {}
            self.flow_in = pyunits.convert(self.flow_vol_in[time], to_units=(pyunits.Mgallons / pyunits.day))
            self.base_fixed_cap_cost = 2.5
            self.mf_cap = self.base_fixed_cap_cost * self.flow_in
            return self.mf_cap
        if self.cost_method == 'wtrnet':
            self.chem_dict = {}
            self.flow_in = pyunits.convert(self.flow_vol_in[time], to_units=(pyunits.m**3 / pyunits.day))
            self.base_fixed_cap_cost = 5.764633
            self.cap_cost_exp = 0.6
            self.mf_cap = (self.base_fixed_cap_cost * self.flow_in ** self.cap_cost_exp) * 1E-3
            self.costing.other_var_cost = (0.015008 * self.flow_in ** 1.072667) * 1E-3
            return self.mf_cap


    def elect(self):
        self.electricity_intensity = Var(initialize=0.18,
                                    bounds=(0, None),
                                    units=pyunits.kWh/pyunits.m**3,
                                    doc='MF electricity intensity [kWh/yr]')
        if self.cost_method == 'twb':
            self.electricity_intensity.fix(0.18)
            return self.electricity_intensity
        if self.cost_method == 'wtrnet':
            self.electricity_intensity_constr = \
                    Constraint(expr=self.electricity_intensity == (91.28175 * self.flow_in ** 0.999957) / 
                                pyunits.convert(self.flow_in, to_units=pyunits.m**3/pyunits.yr))
            return self.electricity_intensity

    def get_costing(self, unit_params=None, year=None):
        '''
        Initialize the unit in WaterTAP3.
        '''
        
        if 'cost_method' in self.unit_params.keys():
            self.cost_method = self.unit_params['cost_method']
            if self.cost_method not in ['twb', 'wtrnet']:
                self.cost_method = 'twb'
        else:
            self.cost_method = 'twb'
        if self.cost_method == 'twb':
            self.basis_year = 2014
        if self.cost_method == 'wtrnet':
            self.basis_year = 2006
        financials.create_costing_block(self, self.basis_year, tpec_or_tic)

        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                                                           doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.elect(),
                                      doc='Electricity intensity [kwh/m3]')
        financials.get_complete_costing(self.costing)