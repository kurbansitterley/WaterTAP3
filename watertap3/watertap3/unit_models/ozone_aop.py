from pyomo.environ import Var, Constraint, Expression, value, units as pyunits
from watertap3.utils import financials
from watertap3.core.wt3_unit_siso import WT3UnitProcessSISO
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

## REFERENCES
## CAPITAL:
# Regression based on Texas Water Board - IT3PR documentation Table 3.24
# https://www.twdb.texas.gov/publications/reports/contracted_reports/doc/1348321632_manual.pdf
## ELECTRICITY:
# "A Review of Ozone Systems Costs for Municipal Applications. Report by the Municipal Committee – IOA Pan American Group"
# (2018) B. Mundy et al. https://doi.org/10.1080/01919512.2018.1467187
# From source: "5.0–5.5 kWh/lb might be used to include total energy consumption of ozone generator, ozone destruct, nitrogen feed system, cooling water pumps, and other miscellaneous
# energy uses."

module_name = 'ozone_aop'

class UnitProcess(WT3UnitProcessSISO):

    def ozone_aop_setup(self):
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time], 
            to_units=pyunits.Mgallons/pyunits.day)

        self.ox_dose = Var(initialize=1,
            units=pyunits.mg/pyunits.L,
            bounds=(0, None),
            doc='Oxidant dose [mg/L]')

        if 'toc_method' in self.unit_params.keys():
            self.toc_method = self.unit_params['toc_method']
            if self.toc_method not in ['user', 'source', 'deprecated']:
                self.toc_method = 'source'
        else:
            self.toc_method = 'source'
        if self.toc_method == 'source':
            try:
                self.toc_in = (self.parent_block().source_df.loc['toc'].value * 1000) * \
                    (pyunits.mg/pyunits.liter)
            except KeyError:
                print('TOC not found in source file. Assuming TOC = 5 mg/L')
                self.toc_in = 5 * (pyunits.mg/pyunits.liter)
        if self.toc_method == 'user':
            self.toc_in = self.unit_params['toc_in'] * (pyunits.mg/pyunits.liter)
        if self.toc_method == 'deprecated':
            self.toc_in = pyunits.convert(self.conc_mass_in[0, 'toc'],
                to_units=pyunits.mg/pyunits.liter)
        if 'aop' in self.unit_params.keys():
            self.aop = self.unit_params['aop']
        else:
            self.aop = False
        if 'contact_time' in self.unit_params.keys():
            self.contact_time = self.unit_params['contact_time'] * pyunits.minutes
        else:
            self.contact_time = 5 * pyunits.minutes
        if 'ct' in self.unit_params.keys():
            self.ct = self.unit_params['ct'] * ((pyunits.mg * pyunits.minute) / (pyunits.liter))
        else:
            self.ct = 1 * ((pyunits.mg * pyunits.minute) / (pyunits.liter))
        if 'mass_transfer' in self.unit_params.keys():
            self.mass_transfer = self.unit_params['mass_transfer'] * pyunits.dimensionless
        else:
            self.mass_transfer = 0.8 * pyunits.dimensionless
        
        self.ozone_consumption = (self.toc_in + self.ct / self.contact_time) / self.mass_transfer
        self.o3_toc_ratio = 1 + ((self.ct / self.contact_time) / self.toc_in)

        if self.aop:
            if 'dose' in self.unit_params.keys():
                self.ox_dose.fix(value(pyunits.convert(self.unit_params['dose']*pyunits.mg/pyunits.L,
                    to_units=pyunits.kg/pyunits.m**3)))
            elif 'dose' not in self.unit_params.keys():
                self.ox_dose_constr = Constraint(
                    expr=self.ox_dose == pyunits.convert((0.5 * self.o3_toc_ratio * self.toc_in), 
                    to_units=(pyunits.kg/pyunits.m**3)))
            else:
                self.ox_dose = 0
            if 'chemical_name' in self.unit_params.keys():
                self.chem_name = self.unit_params['chemical_name']
            else:
                self.chem_name = 'Hydrogen_Peroxide'
            self.chem_dict = {self.chem_name: self.ox_dose}
            self.h2o2_base_cap = 1228
            self.h2o2_cap_exp = 0.2277
        else:
            self.ox_dose.fix(0)
            

    def fixed_cap(self):
        '''
        Fixed capital for Ozone/Ozone AOP unit.

        '''
        self.a, self.b = self.interp_cost_at_dose(self.ozone_consumption())
        self.ozone_cap = self.a * self.flow_in ** self.b
        if self.aop:
            h2o2_flow = self.solution_vol_flow()
            h2o2_cap = self.h2o2_base_cap * h2o2_flow ** self.h2o2_cap_exp
        else:
            h2o2_cap = 0
        ozone_aop_cap = (self.ozone_cap + h2o2_cap) * 1E-3

        return ozone_aop_cap

    def interp_cost_at_dose(self, dose):
        '''
        Determine a, b costing parameters as a function of ozone dose
        :param flow: Volumetric flow into unit
        :param dose: ozone dose
        :return: capital cost
        '''

        def power(flow, a, b):
            return a * flow ** b

        df = pd.read_csv('data/ozone_cost.csv', header=0)

        doses = [1, 5, 10, 15, 20, 25]
        self.flow_interp = []
        self.cost_interp = []
        for i, k in enumerate(df['flow (mgd)']):
            costs = []
            for d in doses:
                cutal_ox = df[str(d)].to_numpy()
                costs.append(cutal_ox[i])
            self.flow_interp.append(k)
            self.cost_interp.append(np.interp(dose, doses, costs))
        (a, b), _ = curve_fit(power, self.flow_interp, self.cost_interp, \
            bounds=[[1E-5, 1E-5], [100000, 5]])
        return a, b

    def elect(self):
        '''
        Electricity intensity for Ozone/Ozone AOP unit.

        '''

        if 'EEO' in self.unit_params.keys():
            self.EEO = Var(initialize=0.1,
                units=(pyunits.kW*pyunits.hr)/(pyunits.m**3),
                bounds=(0, None),
                doc='Electric energy per order of removal (EE/O) [kWh/m3]')
            self.EEO.fix(self.unit_params['EEO'])
            return self.EEO
        else:
            self.ozone_power = Var(initialize=0.1,
                units=pyunits.kW,
                bounds=(0, None),
                doc='Ozone power [kW]')
            self.ozone_power_coeff = Var(initialize=5,
                units=pyunits.kW/(pyunits.lb/pyunits.hr),
                bounds=(0, None),
                doc='Ozone power coeff [kW/lb/hr]')
            self.ozone_power_coeff.fix(5)

            self.ozone_flow = Var(initialize=1,
                    units=pyunits.lb/pyunits.hr,
                    bounds=(0, None),
                    doc='Ozone flow [lb/hr]')
            self.ozone_flow_constr = Constraint(expr=
                    self.ozone_flow == pyunits.convert(self.flow_in * self.ozone_consumption,
                        to_units=pyunits.lb/pyunits.hr))
            self.ozone_power_constr = Constraint(expr=
                            self.ozone_power == self.ozone_power_coeff * self.ozone_flow)
            flow_in_m3hr = pyunits.convert(self.flow_in, 
                                to_units=(pyunits.m**3/pyunits.hour))
            electricity = self.ozone_power / flow_in_m3hr
            return electricity

    def solution_vol_flow(self):
        '''
        Determine oxidant solution flow rate [lb/day]

        :return: Oxidant solution flow [lb/day]
        '''
        chemical_rate = self.flow_in * self.ox_dose  
        chemical_rate = pyunits.convert(chemical_rate,
            to_units=(pyunits.lb/pyunits.day))
        return chemical_rate

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2014
        self.ozone_aop_setup()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.elect(),
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year)