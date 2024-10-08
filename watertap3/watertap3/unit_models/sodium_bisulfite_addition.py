from pyomo.environ import Var, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.core.wt3_unit_pt import WT3UnitProcessPT

## REFERENCE
## CAPITAL:
# Based on costs for SULFURIC ACID ADDITION 93% SOLUTION - FIGURE 5.5.11
# Cost Estimating Manual for Water Treatment Facilities (McGivney/Kawamura) (2008)
# DOI:10.1002/9780470260036
## ELECTRICITY:

module_name = 'sodium_bisulfite_addition'

class UnitProcess(WT3UnitProcessPT):

    def fixed_cap(self):
        '''
        **"unit_params" are the unit parameters passed to the model from the input sheet as a Python dictionary.**

        **EXAMPLE: {'dose': 10}**

        Fixed capital for sodium bisulfite addition is a function of sodium bisulfite dose, sodium bisulfite solution flow, and the number of units.

        :param dose: Sodium bisulfite dose [mg/L]
        :type dose: float
        :return: Sodium bisulfite addition fixed capital cost [$MM]
        '''
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.hr)
        self.number_of_units = 2
        self.base_fixed_cap_cost = 900.97
        self.cap_scaling_exp = 0.6179
        chem_name = 'Sodium_Bisullfite_NaHSO3'
        self.dose = Var(initialize=1,
            bounds=(0, None),
            units=pyunits.kg/pyunits.m**3,
            doc='Dose [kg/m3]')
        self.dose.fix(0.001)
        if 'dose' in self.unit_params.keys():
            self.dose.fix(self.unit_params['dose'] * 1E-3)
        self.chem_dict = {chem_name: self.dose}
        source_cost = self.base_fixed_cap_cost * self.solution_vol_flow() ** self.cap_scaling_exp
        bisulfite_cap = (source_cost * self.tpec_tic * self.number_of_units) * 1E-6
        return bisulfite_cap

    def elect(self):
        '''
        Electricity intensity.

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
        Chemical solution flow in gal/day

        :param solution_density: Solution density [kg/m3]
        :type solution_density: float

        :return: Sodium bisulfite solution flow [gal/day]
        '''
        self.solution_density = 1480 * (pyunits.kg/pyunits.m**3)
        chemical_rate = self.flow_in * self.dose
        chemical_rate = pyunits.convert(chemical_rate,
            to_units=(pyunits.kg/pyunits.day))
        soln_vol_flow = chemical_rate / self.solution_density
        soln_vol_flow = pyunits.convert(soln_vol_flow,
            to_units=pyunits.gallon/pyunits.day)
        return soln_vol_flow  # gal/day

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