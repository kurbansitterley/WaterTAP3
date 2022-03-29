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
tpec_or_tic = 'TPEC'

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
        self.chem_dict = {}
        if self.cost_method == 'twb':
            self.flow_in = pyunits.convert(self.flow_vol_in[time], 
                to_units=(pyunits.Mgallons / pyunits.day))
            self.nf_mem_equipment.fix(0.75)
            self.nf_equip_multiplier.fix(5)
            self.nf_cap_exp.fix(1)
            self.nf_cap_base_constr = Constraint(expr=self.nf_cap_base == 
                        self.nf_mem_equipment * self.nf_equip_multiplier)
        if self.cost_method == 'wtrnet':
            self.flow_in = pyunits.convert(self.flow_vol_in[time], 
                to_units=(pyunits.m**3 / pyunits.day))
            self.nf_cap_base.fix(1.012361 * 1E-3)
            self.nf_cap_exp.fix(0.844997)
            self.costing.other_var_cost = (0.001879 * self.flow_in ** 1.353971) * 1E-3
        self.nf_cap_constr = Constraint(expr=self.nf_fixed_cap == 
                    self.nf_cap_base * self.flow_in ** self.nf_cap_exp)
        return self.nf_fixed_cap

    def elect(self):
        self.electricity_intensity = Var(
                    initialize=0.18,
                    bounds=(0, None),
                    units=pyunits.kWh/pyunits.m**3,
                    doc='MF electricity intensity [kWh/yr]')
        if self.cost_method == 'twb':
            self.electricity_intensity.fix(0.18)
            return self.electricity_intensity
        if self.cost_method == 'wtrnet':
            self.electricity_intensity_constr = \
                    Constraint(expr=self.electricity_intensity ==
                    (164.2818 * self.flow_in ** 0.999976) / 
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
            self.water_recovery.fix(0.83)
            self.basis_year = 2006
        financials.create_costing_block(self, self.basis_year, tpec_or_tic)

        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                                                           doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.elect(),
                                      doc='Electricity intensity [kwh/m3]')
        financials.get_complete_costing(self.costing)