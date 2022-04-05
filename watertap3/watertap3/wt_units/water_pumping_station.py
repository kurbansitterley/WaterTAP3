from pyomo.environ import Expression, value, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE: McGivney & Kawamura (2008) - Figure 5.5.35b & Figure 5.5.36b

module_name = 'water_pumping_station'

class UnitProcess(WT3UnitProcess):

    def fixed_cap(self):
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=(pyunits.Mgallons / pyunits.day))
        if 'pump_type' in self.unit_params.keys():
            self.pump_type = self.unit_params['pump_type']
            if self.pump_type not in ['raw', 'treated']:
                self.pump_type = 'treated'
        else:
            self.pump_type = 'treated'
        if self.pump_type == 'raw':
            self.a = 19370.357574406607
            self.b = 0.9148641590272578
        if self.pump_type == 'treated':
            self.a = 40073.42661387725
            self.b = 0.866701037568153
        pumping_cap = self.tpec_tic * self.a * self.flow_in ** self.b * 1E-6
        return pumping_cap

    def elect(self):
        self.pump_eff = 0.9
        self.motor_eff = 0.9

        flow_in_m3hr = pyunits.convert(self.flow_in,
            to_units=(pyunits.m**3/pyunits.hr))
        flow_in_gpm = value(pyunits.convert(self.flow_in,
            to_units=(pyunits.gallons/pyunits.minute)))

        if 'lift_height' in self.unit_params.keys():
            self.lift_height = self.unit_params['lift_height']
        else:
            self.lift_height = 100 # 100 ft = 3 bar of dynamic head, assume 90% efficiency for motor and pump

        if 'pump_power' in self.unit_params.keys():
            self.pump_power_hp = self.unit_params['pump_power'] * pyunits.hp
            self.pump_power_kw = pyunits.convert(self.pump_power_hp,
                to_units=pyunits.kilowatts)
        else:
            self.pump_power_kw = (0.746 * flow_in_gpm * self.lift_height / \
                (3960 * self.pump_eff * self.motor_eff)) * pyunits.kilowatts

        self.elect_intens = self.pump_power_kw / flow_in_m3hr
        return self.elect_intens

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2007
        tpec_tic = 'TPEC'
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.elect(),
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)