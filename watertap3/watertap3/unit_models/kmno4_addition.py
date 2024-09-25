from pyomo.environ import Var, Constraint, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.core.wt3_unit_pt import WT3UnitProcessPT

## REFERENCE
## CAPITAL:
# McGivney/Kawamura (2008) Fig. 5.5.16
## ELECTRICITY:
# citation here

module_name = 'kmno4_addition'

class UnitProcess(WT3UnitProcessPT):

    def kmno4_setup(self):

        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time], 
            to_units=pyunits.m**3/pyunits.hr)

        self.dose = Var(initialize=1,
            bounds=(0, None),
            units=pyunits.kg/pyunits.m**3,
            doc='Dose [kg/m3]')
        self.dose.fix(0.003)
        if 'dose' in self.unit_params.keys():
            self.dose.fix(self.unit_params['dose'] * 1E-3)

        self.kmno4_density = Var(initialize=1,
            bounds=(0, None),
            units=pyunits.kg/pyunits.m**3,
            doc='KMnO4 solution density [kg/m3]')
        self.kmno4_density.fix(1020)  #CAIROX potassium permanganate, carusllc.com

        self.kmno4_ratio_in_soln = Var(initialize=0.03,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc='KMnO4 ratio in solution')
        self.kmno4_ratio_in_soln.fix(0.03)

        self.kmno4_capital_A = Var(initialize=21434,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='KMnO4 capital A factor')
        self.kmno4_capital_A.fix(21434)

        self.kmno4_capital_B = Var(initialize=0.0758,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='KMnO4 capital B factor')
        self.kmno4_capital_B.fix(0.0758)

        self.motor_eff = Var(initialize=0.9,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc='Pump efficiency')
        self.motor_eff.fix(0.9)

        self.pump_eff = Var(initialize=0.9,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc='Motor efficiency')
        self.pump_eff.fix(0.9)

        self.lift_height = Var(initialize=100,
            bounds=(0, None),
            units=pyunits.ft,
            doc='Pump lift height')
        self.lift_height.fix(100)

        self.kmno4_feed_rate = Var(initialize=500,
            bounds=(0, None),
            units=pyunits.lb/pyunits.day,
            doc='KMnO4 feed rate [lb/day]')

        self.kmno4_soln_flow= Var(initialize=500,
            bounds=(0, None),
            units=pyunits.gallon/pyunits.min,
            doc='KMnO4 solution feed rate [gal/min]')

        self.kmno4_feed_rate_constr = Constraint(expr=
            self.kmno4_feed_rate == pyunits.convert(
                self.flow_in * self.dose, 
                to_units=pyunits.lb/pyunits.day
            ))

        self.kmno4_soln_flow_constr = Constraint(expr=
            self.kmno4_soln_flow == pyunits.convert(pyunits.convert(
                self.kmno4_feed_rate, to_units=pyunits.kg/pyunits.day
            ) / self.kmno4_density / self.kmno4_ratio_in_soln, 
            to_units=pyunits.gallon/pyunits.minute)
            )
        
        self.chem_dict = {'Potassium_Permanganate_KMnO4': self.dose}

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2007
        tpec_tic = 'TPEC'
        self.kmno4_setup()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
            (self.kmno4_capital_A * self.kmno4_feed_rate ** self.kmno4_capital_B) *
            self.tpec_tic * 1E-6,
            doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=
            ((0.746 * self.kmno4_soln_flow * self.lift_height) / 
            (3960 * self.pump_eff * self.motor_eff)) / self.flow_in, 
            doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)