from pyomo.environ import Expression, Var, Constraint, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCES: 
# Capital: 
#   'twb': Texas Water Board User's Manual for Integrated Treatment Train Toolbox - Potable Reuse (IT3PR) Version 2.0.
#    'wtrnet': Joksimovic, D. (2006). Decision Support System for Planning of integrated Water Reuse Projects. (PhD Thesis).
# Electricity
#   'twb': Plappally, A. K., & Lienhard V, J. H. (2012). doi:10.1016/j.rser.2012.05.022
#   'wtrnet': Joksimovic, D. (2006). Decision Support System for Planning of integrated Water Reuse Projects. (PhD Thesis).

module_name = 'nanofiltration'

class UnitProcess(WT3UnitProcess):

    def fixed_cap(self):
        time = self.flowsheet().config.time.first()

        self.nf_fixed_cap = Var(initialize=10,
            bounds=(0, None),
            doc='NF Fixed Capital')
        self.nf_cap_base = Var(initialize=1,
            bounds=(0, None),
            doc='NF Capital Basis')
        self.nf_mem_equipment = Var(initialize=0.5,
            bounds=(0, None),
            doc='NF Membrane Equipment [$M/MGD]')
        self.nf_equip_multiplier = Var(initialize=1,
            bounds=(0, None),
            doc='NF Equipment Multiplier [$/MGD]')
        self.nf_cap_exp = Var(initialize=0.1,
            bounds=(0, None),
            doc='NF Capital Exponent')
        
        if self.cost_method == 'twb':
            self.flow_in = pyunits.convert(self.flow_vol_in[time], 
                to_units=(pyunits.Mgallons / pyunits.day))
            self.nf_mem_equipment.fix(0.75)
            self.nf_equip_multiplier.fix(5)
            self.nf_cap_exp.fix(0.7)
            self.nf_cap_base_constr = Constraint(expr=self.nf_cap_base == 
                    self.nf_mem_equipment * self.nf_equip_multiplier)
            self.nf_cap_constr = Constraint(expr=self.nf_fixed_cap == 
                    self.nf_cap_base * self.flow_in ** self.nf_cap_exp)
        if self.cost_method == 'wtrnet':
            self.flow_in = pyunits.convert(self.flow_vol_in[time], 
                to_units=(pyunits.m**3 / pyunits.day))
            self.nf_cap_base.fix(1.012361 * 1E-3)
            self.nf_cap_exp.fix(0.844997)
            self.costing.other_var_cost = (0.001879 * self.flow_in ** 1.353971) * 1E-3
            self.nf_cap_constr = Constraint(expr=self.nf_fixed_cap == self.tpec_tic *
                        (self.nf_cap_base * self.flow_in ** self.nf_cap_exp))

    def elect(self):
        time = self.flowsheet().config.time.first()
        
        self.flow_in = pyunits.convert(self.flow_vol_in[time], 
                    to_units=pyunits.m**3/pyunits.hr)
        
        self.electricity_intensity = Var(
            initialize=0.18,
            bounds=(0, None),
            units=pyunits.kWh/pyunits.m**3,
            doc='NF electricity intensity [kWh/m3]')
        
        self.pressure = Var(
            initialize=4,
            bounds=(0, 20),
            units=pyunits.bar,
            doc='NF operating pressure [bar]')
        self.pressure.fix(10)

        self.pump_eff = Var(
            initialize=0.8,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc='NF pump efficiency ')
        self.pump_eff.fix(0.8)

        for k, v in self.unit_params.items():
            if k in ['pressure', 'pump_eff']:
                getattr(self, k).fix(v)

        self.pump_power = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.kW,
            doc='NF pump power required [kW]')
        
        self.pump_power_constr = Constraint(expr=self.pump_power == 
                (pyunits.convert(self.flow_in * 
                pyunits.convert(self.pressure, to_units=pyunits.Pa),
                to_units=pyunits.kW)) / self.pump_eff)
        
        self.electricity_intensity_constr = Constraint(expr=
            self.electricity_intensity == self.pump_power /
            self.flow_in)

        # if self.cost_method == 'twb':
        #     self.electricity_intensity.fix(0.18)
            # return self.electricity_intensity
        # if self.cost_method == 'wtrnet':
        #     self.electricity_intensity_constr = \
        #             Constraint(expr=self.electricity_intensity ==
        #             (164.2818 * self.flow_in ** 0.999976) / 
        #             pyunits.convert(self.flow_in, to_units=pyunits.m**3/pyunits.yr))
            # return self.electricity_intensity
        # else:

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        tpec_tic = 'TIC'
        if 'cost_method' in self.unit_params.keys():
            self.cost_method = self.unit_params['cost_method']
            if self.cost_method not in ['twb', 'wtrnet']:
                self.cost_method = 'twb'
        else:
            self.cost_method = 'twb'
        if self.cost_method == 'twb':
            self.basis_year = 2014
        if self.cost_method == 'wtrnet':
            # self.water_recovery.fix(0.83)
            self.basis_year = 2006
        
        self.fixed_cap()
        self.elect()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.nf_fixed_cap,
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.electricity_intensity,
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=self.basis_year, tpec_tic=tpec_tic)