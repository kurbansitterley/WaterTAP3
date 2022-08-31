from pyomo.environ import Expression, Var, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit_pt import WT3UnitProcessPT

## REFERENCE: Derived from Voutchkov (2018) Table 4.6 and 4.7

module_name = 'well_field'

class UnitProcess(WT3UnitProcessPT):

    def fixed_cap(self):
        self.pipe_cost_factor_dict = {'emwd': 82600}
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.hr)
        self.base_fixed_cap_cost = 4731.6
        self.cap_scaling_exp = 0.9196
        try:
            self.pipe_distance = self.unit_params['pipe_distance'] * pyunits.miles
            self.pipe_diameter = 8 * pyunits.inches
            try:
                self.pipe_cost_case = self.unit_params['pipe_cost_case']
                self.pipe_cost_basis = self.pipe_cost_factor_dict[self.pipe_cost_case]
            except:
                self.pipe_cost_basis = 35000
            self.pipe_fixed_cap_cost = (self.pipe_cost_basis * \
                self.pipe_distance * self.pipe_diameter)
            well_cap = (self.base_fixed_cap_cost * self.flow_in ** self.cap_scaling_exp + \
                self.pipe_fixed_cap_cost) * 1E-6
            return well_cap
        except:
            well_cap = self.base_fixed_cap_cost * self.flow_in ** self.cap_scaling_exp * 1E-6
            return well_cap

    def elect(self):
        time = self.flowsheet().config.time
        t = time.first()
        try:
            self.pump = self.unit_params['pump']
            if self.pump not in ['yes', 'no']:
                self.pump = 'yes'
        except (KeyError, TypeError) as e:
            self.pump = 'yes'
        self.lift_height = Var(time,
            initialize=100,
            bounds=(0, 1E5),
            units=pyunits.ft,
            doc='Lift height for well pump [ft]')
        self.pump_eff = 0.9 * pyunits.dimensionless
        self.motor_eff = 0.9 * pyunits.dimensionless
        if self.pump == 'yes':
            if 'lift_height' in self.unit_params.keys():
                self.lift_height.fix(self.unit_params['lift_height'])
            else:
                self.lift_height.fix(100)
            flow_in_gpm = pyunits.convert(self.flow_vol_in[t],
                to_units=pyunits.gallons/pyunits.minute)
            electricity = (0.746 * flow_in_gpm * self.lift_height[t] / \
                (3960 * self.pump_eff * self.motor_eff)) / self.flow_in
            return electricity
        else:
            return 0

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2018
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.elect(),
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year)