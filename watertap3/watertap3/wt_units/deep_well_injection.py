from pyomo.environ import Expression, Var, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit_pt import WT3UnitProcessPT

## REFERENCE:
## CAPITAL
# Developed from Kay Bailey and produced water case study data

module_name = 'deep_well_injection'

class UnitProcess(WT3UnitProcessPT):

    def fixed_cap(self):
        '''
        Fixed capital cost for deep well injection.

        :param unit_params: Input parameter dictionary from input sheet.
        :type unit_params: dict
        :param lift_height: Lift height for pump [ft]
        :type lift_height: float
        :param pipe_distance: Piping distance to deep well injection site
        :type pipe_distance: float
        :return: Fixed capital cost for deep well injection [$MM]
        '''
        time = self.flowsheet().config.time
        t = self.flowsheet().config.time.first()
        self.lift_height = Var(time,
            initialize=400,
            units=pyunits.ft,
            doc='Lift height for pump [ft]')
        self.lift_height.fix(400)
        if 'lift_height' in self.unit_params.keys():
            self.lift_height.fix(self.unit_params['lift_height'])
        self.flow_in = pyunits.convert(self.flow_vol_in[t],
            to_units=pyunits.m**3/pyunits.hr)
            
        self.delta_pressure = self.lift_height[t] * 0.0299
        self.cap_scaling_exp = 0.7
        self.cap_scaling_val = 473.2
        self.well_pump_fixed_cap_cost = 16.9  # this is wells/pumps fixed capital AFTER applying TIC factor -- DOES NOT INCLUDE ANY PIPING
        self.pipe_cost_basis = 35000  # $ / (inch * mile) -- this value taken from produced water case studies in WT3 Excel model

        if 'pipe_distance' in self.unit_params.keys():
            self.pipe_distance = self.unit_params['pipe_distance'] * pyunits.miles
        else:
            self.pipe_distance = 0 * pyunits.miles
        self.pipe_diameter = 8 * pyunits.inches
        
        self.pipe_fixed_cap_cost = (self.pipe_cost_basis * self.pipe_distance * self.pipe_diameter) * 1E-6
        self.tot_fixed_cap = self.well_pump_fixed_cap_cost + self.pipe_fixed_cap_cost
        cap_scaling_factor = self.flow_in / self.cap_scaling_val
        deep_well_cap = self.tot_fixed_cap * cap_scaling_factor ** self.cap_scaling_exp
        return deep_well_cap

    def elect(self):
        '''
        Electricity intensity for deep well injection [kWh/m3]

        :param lift_height: Lift height for pump [ft]
        :type lift_height: float
        :return: Electricity intensity [kWh/m3]
        '''
        t = self.flowsheet().config.time.first()
        self.pump_eff = 0.9
        self.motor_eff = 0.9
        flow_in_gpm = pyunits.convert(self.flow_in,
            to_units=(pyunits.gallon/pyunits.minute))
        electricity = (0.746 * flow_in_gpm * self.lift_height[t] / \
            (3960 * self.pump_eff * self.motor_eff)) / self.flow_in
        return electricity

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2011
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.elect(),
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year)