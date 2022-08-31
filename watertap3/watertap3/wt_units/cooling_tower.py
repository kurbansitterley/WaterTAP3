from pyomo.environ import Expression, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE: ADD REFERENCE HERE

module_name = 'cooling_tower'

class UnitProcess(WT3UnitProcess):

    def fixed_cap(self):
        time = self.flowsheet().config.time.first()
        flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.hr)
        atmospheric_pressure = 101325  # Pa
        ambient_temperature = 20
        relative_humidity = 0.5
        self.cp = 4.18  # heat capacity
        self.cycles = 5.0
        self.ci_circ = 0.0  # unset for the ratio, but could be defined in future
        self.ci_makeup = 0.0  # unset for the ratio, but could be defined in future
        self.latent_heat = 2264.76  # latent heat of vaporization MJ/m3
        self.evap_fraction = 0.85  # fraction of total heat rejected by latent heat transfer - 0.9 in EPRI report
        self.approach = 5.55  # see Miara & Vorosmarty 2017 for assumption
        self.wet_bulb_temp = 20  # unless specified
        self.ttd = 4.44  # Celsius based on typical range of 8F from EPRI 2004
        self.range = 11.11  # Celsius based on typical range of 20F from EPRI 2004
        if 'cycles' in self.unit_params.keys():
            # print('Cycles of concentration are set based on parameters to:',self.unit_params['cycles'])
            self.cycles = self.unit_params['cycles'];
        else:
            self.cycles = 5.0;
            print('If cycles are used, assuming cycles of concentration is:', self.cycles)
        ## EVPORATION FRACTION PROVIDED BY USER
        if self.unit_params['method'] == 'evaporation_fraction':
            self.make_up = self.flow_vol_in[time]
            self.evaporation = self.unit_params['evaporation_fraction'] * self.flow_vol_in[time]
            self.blowdown = self.make_up - self.evaporation
        ## MASS BALANCE APPROACH
        if self.unit_params['method'] == 'make_up_mass_balance':
            self.make_up = self.flow_vol_in[time]
            self.blowdown = self.make_up / self.cycles
            self.evaporation = self.make_up - self.blowdown  # evaporation assumed to go to waste outlet (out of system) should not go to surface discharge
        ## MASS BALANCE AND PLANT INFO USED TO CALCULATE CYCLES OF CONCENTRATION, BLOWDOWN, EVAPORATION
        if self.unit_params['method'] == 'plant_info_with_makeup':
            print('Make up is given, cycles of concentration are calculated as a result of assumptions')
            self.set_plant_info(self)
            self.nameplate = self.unit_params['nameplate']
            self.heat_in = self.nameplate / self.eff
            self.desired_heat_cond = self.heat_in - (self.heat_in * self.eff) - \
                (self.heat_in * self.heat_sink)
            self.evaporation = self.desired_heat_cond * (self.evap_fraction / self.latent_heat)
            self.flow_vol_in.unfix()
            self.unfix_inlet_to_train(self)
            self.flow_vol_waste.fix(self.evaporation)
            self.make_up = self.unit_params['make_up']
            self.blowdown = self.make_up - self.evaporation
            self.cycles = self.make_up / self.blowdown
        ## MASS BALANCE AND PLANT INFO USED TO CALCULATE CYCLES OF CONCENTRATION, BLOWDOWN, EVAPORATION
        if self.unit_params['method'] == 'plant_info_without_makeup':
            print('make up not given as a result of assumptions, cycles of concentration are given')
            self.set_plant_info(self)
            self.nameplate = self.unit_params['nameplate']
            self.heat_in = self.nameplate / self.eff
            self.desired_heat_cond = self.heat_in - (self.heat_in * self.eff) - \
                (self.heat_in * self.heat_sink)
            self.evaporation = self.desired_heat_cond * (self.evap_fraction / self.latent_heat)
            self.flow_vol_in.unfix()
            self.unfix_inlet_to_train(self)
            self.flow_vol_waste.fix(self.evaporation)
            self.make_up = self.flow_vol_in[time]
            self.blowdown = self.make_up - self.evaporation
            self.make_up = self.flow_vol_in[time]
            self.blowdown = self.make_up / self.cycles
            self.evaporation = self.make_up - self.blowdown  # evaporation assumed to go to waste outlet (out of system) should not go to surface discharge
        self.water_recovery.fix(self.blowdown / self.make_up)
        
        ct_cap = self.flow_vol_in[time] * 1E-9  # $M
        return ct_cap

    def elect(self):
        electricity = 1E-9
        return electricity

    def unfix_inlet_to_train(self):
        if hasattr(self.parent_block(), 'pfd_dict'):
            for key in self.parent_block().pfd_dict:
                if self.parent_block().pfd_dict[key]['type'] == 'intake':
                    print('unfixing intake:', key)
                    getattr(self.parent_block(), key).flow_vol_in.unfix()
        else:
            print('assuming test with source1')
            self.parent_block().source1.flow_vol_in.unfix()

    def set_plant_info(self):
        if self.unit_params['fuel'] == 'nuclear':
            self.heat_sink = 0.0
            self.eff = 0.3
        if self.unit_params['fuel'] == 'natural_gas_cc':
            self.heat_sink = 0.2
            self.eff = 0.6
        if self.unit_params['fuel'] == 'coal':
            self.heat_sink = 0.12
            self.eff = 0.35

        if 'efficiency' in self.unit_params.keys():
            print('Thermal efficiency of plant set based on parameters to:', \
                self.unit_params['efficiency'])
            self.eff = self.unit_params['efficiency'];
        else:
            print('Assuming thermal efficiency of plant is:', self.eff)

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2020
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')  # $M
        self.electricity = Expression(expr=self.elect(),
                doc='Electricity intensity [kWh/m3]')  # kwh/m3
        financials.get_complete_costing(self.costing, basis_year=basis_year)