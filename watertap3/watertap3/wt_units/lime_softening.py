from pyomo.environ import Var, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# Minnesota Rural Water Association, Chapter 16 Lime Softening
# #(https://www.mrwa.com/WaterWorksMnl/Chapter%2016%20Lime%20Softening.pdf)
# https://www.necoindustrialwater.com/analysis-ion-exchange-vs-lime-softening/

module_name = 'lime_softening'

class UnitProcess(WT3UnitProcess):

    def fixed_cap(self):
        '''
        **"unit_params" are the unit parameters passed to the model from the input sheet as a Python dictionary.**

        **EXAMPLE: {'lime': 10}**

        Fixed capital for lime addition is a function of lime dose, lime solution flow, and the number of units.

        :param lime: Lime dose [mg/L]
        :type lime: float

        :return: Lime addition fixed capital cost [$MM]
        '''
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.hr)
        self.base_fixed_cap_cost = 0.0704
        self.cap_scaling_exp = 0.7306
        chem_name = 'Lime_Suspension_CaOH_2'
        self.dose = Var(initialize=1,
            bounds=(0, None),
            units=pyunits.kg/pyunits.m**3,
            doc='Dose [kg/m3]')
        self.dose.fix(0.0023)
        if 'dose' in self.unit_params.keys():
            self.dose.fix(self.unit_params['dose'] * 1E-3)
        # self.magnesium_dissolved_lime = pyunits.convert(self.unit_params['lime'] * \
        #     (pyunits.mg/pyunits.liter), 
        #     to_units=(pyunits.kg/pyunits.m**3))
        self.magnesium_dissolved_factor = 30 * pyunits.dimensionless  #
        self.chemical_dosage = self.magnesium_dissolved_factor * self.dose
        self.chem_dict = {chem_name: self.chemical_dosage}
        lime_cap = self.base_fixed_cap_cost * self.flow_in ** self.cap_scaling_exp
        return lime_cap

    def elect(self):
        '''
        Electricity intensity for chemical additions is a function of lift height, pump efficiency, and motor efficiency.


        :return: Electricity intensity [kWh/m3]
        '''
        self.lift_height = 100 * pyunits.ft
        self.pump_eff = 0.9 * pyunits.dimensionless
        self.motor_eff = 0.9 * pyunits.dimensionless
        soln_vol_flow = pyunits.convert(self.solution_vol_flow(),
            to_units=(pyunits.gallon/pyunits.minute))
        electricity = (0.746 * soln_vol_flow * self.lift_height / \
            (3960 * self.pump_eff * self.motor_eff)) / self.flow_in
        return electricity

    def solution_vol_flow(self):
        '''
        Determine alum solution flow rate in gal / day


        :return: Lime solution flow [gal/day]
        '''
        self.solution_density = 1250 * (pyunits.kg/pyunits.m**3)
        chemical_rate = self.flow_in * self.chemical_dosage
        chemical_rate = pyunits.convert(chemical_rate,
            to_units=(pyunits.kg/pyunits.day))
        soln_vol_flow = chemical_rate / self.solution_density
        return soln_vol_flow

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2020
        tpec_tic = 'TPEC'
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.elect(),
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)