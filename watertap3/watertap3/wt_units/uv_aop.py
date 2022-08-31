import pandas as pd
from pyomo.environ import Expression, Constraint, Var, Param, units as pyunits
from scipy.optimize import curve_fit
from watertap3.utils import financials
from watertap3.wt_units.wt_unit_siso import WT3UnitProcessSISO

## REFERENCE:
## CAPITAL
# Interpolated values for UV cost based on Texas Water Board - IT3PR documentation Table 3.22
# https://www.twdb.texas.gov/publications/reports/contracted_reports/doc/1348321632_manual.pdf
## ELECTRICITY

module_name = 'uv_aop'

class UnitProcess(WT3UnitProcessSISO):

    def fixed_cap(self):
        '''
        **"unit_params" are the unit parameters passed to the model from the input sheet as a Python dictionary.**

        **EXAMPLE: {'aop': True, 'uv_dose': 350, 'dose': 5, 'chemical_name': 'Hydrogen Peroxide'}**

        :param aop: (**required**) Boolean that determines if UV is used with AOP. Must be either True or False
        :type aop: bool
        :param uvt_in: (**optional**, default is 0.9) UV transmission (UVT) into unit
        :type uvt_in: float
        :param uv_dose: (**optional**, default is 100) Reduction Equivalent Dose (RED)  [mJ/cm2]
        :type uv_dose: float
        :param dose: (**optional**, no default) Dose for oxidant (if AOP) [mg/L]
        :type dose: float
        :param chemical_name: (**optional**, default is ``'Hydrogen_Peroxide'``) Name of oxidant used for AOP.
        :type chemical_name: str


        :return: Fixed capital for UV or UV+AOP unit [$MM]
        '''
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.hr)
        try:
            self.aop = self.unit_params['aop']
        except:
            self.aop = False
        if self.aop:
            try:
                self.dose = self.unit_params['dose'] * (pyunits.mg/pyunits.liter)
            except:
                self.dose = 5 * (pyunits.mg/pyunits.liter)
            try:
                self.chem_name = self.unit_params['chemical_name']
            except:
                self.chem_name = 'Hydrogen_Peroxide'
            self.ox_dose = pyunits.convert(self.dose,
                to_units=(pyunits.kg/pyunits.m**3))
            self.chem_dict = {'Hydrogen_Peroxide': self.ox_dose}
            self.h2o2_base_cap = 1228
            self.h2o2_cap_exp = 0.2277
            self.h2o2_cap = (self.h2o2_base_cap * self.solution_vol_flow() ** self.h2o2_cap_exp) * 1E-3
        else:
            
            self.h2o2_cap = 0
        try:
            self.uvt_in = self.unit_params['uvt_in']
            self.uv_dose = self.unit_params['uv_dose']
        except:
            self.uvt_in = 0.9
            self.uv_dose = 100
        self.a, self.b = self.uv_regress()
        flow_in_mgd = pyunits.convert(self.flow_in,
            to_units=(pyunits.Mgallons/pyunits.day))
        self.uv_cap = (self.a * flow_in_mgd ** self.b) * 1E-3
        self.uv_aop_cap = self.uv_cap + self.h2o2_cap
        return self.uv_aop_cap

    def elect(self):
        
        self.EEO = Var(initialize=0.1,
            units=(pyunits.kW*pyunits.hr)/(pyunits.m**3),
            bounds=(0, 250),
            doc='Electric energy per order of removal (EE/O) [kWh/m3]')
        self.EEO.fix(0.1)
        
        self.lamp_eff = Var(initialize=0.30,
            units=pyunits.dimensionless,
            doc='UV lamp efficiency [dimensionless]')
        self.lamp_eff.fix(0.3)

        self.order_of_mag_removal = Var(initialize=1,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc='Order-of-magnitude removal for contaminant of interest')
        self.order_of_mag_removal.fix(1)

        self.lamp_power = Var(initialize=40,
            units=pyunits.kW,
            doc='UV system lamp power [kW]')
        
        for k, v in self.unit_params.items():
            if k in ['EEO', 'lamp_eff', 'order_of_mag_removal']:
                getattr(self, k).fix(v)

        self.lamp_power_constr = Constraint(expr=self.lamp_power == 
                (self.EEO * self.flow_in * self.order_of_mag_removal) / self.lamp_eff)

        if self.aop:
            self.lift_height = Var(initialize=1,
                units=pyunits.ft,
                bounds=(0, None),
                doc='Lift height for H2O2 pump [ft]')
                            
            self.lift_height.fix(100)

            self.pump_eff = Var(initialize=0.9,
                units=pyunits.dimensionless,
                bounds=(0, 1),
                doc='Pump efficiency for H2O2 pump')

            self.pump_eff.fix(0.9)

            self.motor_eff = Var(initialize=0.9, 
                units=pyunits.dimensionless,
                bounds=(0, 1), 
                doc='Motor efficiency for H2O2 pump')

            self.motor_eff.fix(0.9)

            self.h2o2_density = Param(initialize=1130,
                units=pyunits.kg/pyunits.m**3)

            self.ox_pump_power = Var(initialize=100, 
                units=pyunits.kW,
                bounds=(0, None), 
                doc='Power for H2O2 pump')

            for k, v in self.unit_params.items():
                if k in ['lift_height', 'pump_eff', 'motor_eff']:
                    getattr(self, k).fix(v)
            
            soln_vol_flow = pyunits.convert(self.solution_vol_flow() / self.h2o2_density, 
                                to_units=(pyunits.gallon/pyunits.minute))
            self.ox_power_constr = Constraint(expr=self.ox_pump_power == 
                            (0.746 * soln_vol_flow * self.lift_height / \
                                (3960 * self.pump_eff * self.motor_eff)))
            self.ox_elect = self.ox_pump_power / self.flow_in
        else:
            self.ox_elect = 0

        return (self.EEO * self.order_of_mag_removal) + self.ox_elect

    def uv_regress(self):
        '''
        Determine a, b costing parameters as a function of flow, UVT, and UV dose for unit.

        :return: a, b
        '''

        def power_curve(x, a, b):
            '''
            Return power curve. Used for fitting cost data to determine a, b.
            '''
            return a * x ** b

        self.df = pd.read_csv('data/uv_cost.csv', index_col='flow')
        self.flow_points = [1E-8]
        self.flow_list = [1E-8, 1, 3, 5, 10, 25]  # flow in MGD
        for flow in self.flow_list[1:]:
            temp = self.df.loc[flow]
            cost = temp[((temp.dose == self.uv_dose) & (temp.uvt == self.uvt_in))]
            cost = cost.iloc[0]['cost']
            self.flow_points.append(cost)
        (self.a, self.b), _ = curve_fit(power_curve, self.flow_list, self.flow_points)

        return self.a, self.b

    def solution_vol_flow(self):  # m3/hr
        '''
        Determine oxidant solution flow rate in gal / day

        :return: Oxidant solution flow [gal/day]
        '''
        chemical_rate = self.flow_in * self.ox_dose  # kg/hr
        chemical_rate = pyunits.convert(chemical_rate,
            to_units=(pyunits.lb/pyunits.day))
        soln_vol_flow = chemical_rate
        return soln_vol_flow  # lb / day

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2014
        tpec_tic = 'TPEC'
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')  # $M
        self.electricity = Expression(expr=self.elect(),
                doc='Electricity intensity [kWh/m3]')  # kwh/m3
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)