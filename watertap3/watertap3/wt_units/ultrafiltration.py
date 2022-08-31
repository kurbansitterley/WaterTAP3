from pyomo.environ import Expression, Var, Constraint, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCES: 
# Capital: 
#   'twb': Texas Water Board User's Manual for Integrated Treatment Train Toolbox - Potable Reuse (IT3PR) Version 2.0.
#    'poseidon':
# Electricity
#   'twb': Plappally, A. K., & Lienhard V, J. H. (2012). doi:10.1016/j.rser.2012.05.022
#   'poseidon': Joksimovic, D. (2006). Decision Support System for Planning of integrated Water Reuse Projects. (PhD Thesis).

module_name = 'ultrafiltration'

class UnitProcess(WT3UnitProcess):

    def fixed_cap(self):
        time = self.flowsheet().config.time.first()

        self.uf_fixed_cap = Var(initialize=10,
            bounds=(0, None),
            doc='UF Fixed Capital')
        self.uf_cap_base = Var(initialize=1,
            bounds=(0, None),
            doc='UF Capital Basis')
        self.uf_mem_equipment = Var(initialize=0.5,
            bounds=(0, None),
            doc='UF Membrane Equipment [$M/MGD]')
        self.uf_equip_multiplier = Var(initialize=1,
            bounds=(0, None),
            doc='UF Equipment Multiplier [$/MGD]')
        self.uf_cap_exp = Var(initialize=0.1,
            bounds=(0, None),
            doc='UF Capital Exponent')
        
        if self.cost_method == 'twb':
            self.flow_in = pyunits.convert(self.flow_vol_in[time], 
                to_units=(pyunits.Mgallons / pyunits.day))
            self.uf_mem_equipment.fix(0.5)
            self.uf_equip_multiplier.fix(5)
            self.uf_cap_exp.fix(0.7)
            self.uf_cap_base_constr = Constraint(expr=self.uf_cap_base == 
                self.uf_mem_equipment * self.uf_equip_multiplier)
            self.uf_cap_constr = Constraint(expr=self.uf_fixed_cap == 
                self.uf_cap_base * self.flow_in ** self.uf_cap_exp)
        if self.cost_method == 'poseidon':
            self.flow_in = pyunits.convert(self.flow_vol_in[time], 
                to_units=(pyunits.m**3 / pyunits.day))
            self.uf_cap_base.fix(5.764633 * 1E-3)
            self.uf_cap_exp.fix(0.6)
            # self.costing.other_var_cost = (0.014016 * self.flow_in ** 1.072667) * 1E-3
            self.costing.other_var_cost = (0.014016 * self.flow_in ** 0.7) * 1E-3
            self.uf_cap_constr = Constraint(expr=self.uf_fixed_cap == self.tpec_tic * 
                (self.uf_cap_base * self.flow_in ** self.uf_cap_exp))


    def elect(self):
        time = self.flowsheet().config.time.first()
        
        self.flow_in = pyunits.convert(self.flow_vol_in[time], 
                    to_units=pyunits.m**3/pyunits.hr)
        
        self.electricity_intensity = Var(
            initialize=0.18,
            bounds=(0, None),
            units=pyunits.kWh/pyunits.m**3,
            doc='UF electricity intensity [kWh/m3]')
        
        self.pressure = Var(
            initialize=4,
            bounds=(0, 60),
            units=pyunits.bar,
            doc='UF operating pressure [bar]')
        self.pressure.fix(8)

        self.pump_eff = Var(
            initialize=0.8,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc='UF pump efficiency')
        self.pump_eff.fix(0.8)

        for k, v in self.unit_params.items():
            if k in ['pressure', 'pump_eff']:
                getattr(self, k).fix(v)

        self.pump_power = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.kW,
            doc='UF pump power required [kW]')
        
        self.pump_power_constr = Constraint(expr=self.pump_power == 
                (pyunits.convert(self.flow_in * 
                pyunits.convert(self.pressure, to_units=pyunits.Pa),
                to_units=pyunits.kW)) / self.pump_eff)
        
        self.electricity_intensity_constr = Constraint(expr=
            self.electricity_intensity == self.pump_power /
            self.flow_in)
        # self.electricity_intensity = Var(
        #             initialize=0.18,
        #             bounds=(0, None),
        #             units=pyunits.kWh/pyunits.m**3,
        #             doc='UF electricity intensity [kWh/yr]')
        # if self.cost_method == 'twb':
        #     self.electricity_intensity.fix(0.231344952)
        #     return self.electricity_intensity
        # if self.cost_method == 'poseidon':
        #     self.electricity_intensity_constr = \
        #             Constraint(expr=self.electricity_intensity == (109.5 * self.flow_in) / 
        #             pyunits.convert(self.flow_in, to_units=pyunits.m**3/pyunits.yr))
        #     return self.electricity_intensity

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        tpec_tic = 'TIC'
        if 'cost_method' in self.unit_params.keys():
            self.cost_method = self.unit_params['cost_method']
            if self.cost_method not in ['twb', 'poseidon']:
                self.cost_method = 'poseidon'
        else:
            self.cost_method = 'poseidon'
        if self.cost_method == 'twb':
            self.basis_year = 2014
        if self.cost_method == 'poseidon':
            # self.water_recovery.fix(0.85)
            self.basis_year = 2006
        
        if 'water_recovery' in self.unit_params.keys():
            self.water_recovery.fix(self.unit_params['water_recovery'])

        self.fixed_cap()
        self.elect()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.uf_fixed_cap,
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.electricity_intensity,
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=self.basis_year, tpec_tic=tpec_tic)