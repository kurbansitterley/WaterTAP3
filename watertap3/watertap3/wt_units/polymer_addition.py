from pyomo.environ import Var, Constraint, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# McGivney/Kawamura (2008) Fig. 5.5.16
## ELECTRICITY:
# citation here

module_name = 'polymer_addition'

class UnitProcess(WT3UnitProcess):

    def poly_setup(self):

        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time], 
            to_units=pyunits.m**3/pyunits.hr)

        self.dose = Var(initialize=1,
            bounds=(0, None),
            units=pyunits.kg/pyunits.m**3,
            doc='Dose [kg/m3]')
        self.dose.fix(0.00001)
        if 'dose' in self.unit_params.keys():
            self.dose.fix(self.unit_params['dose'] * 1E-3)

        self.polymer_density = Var(initialize=1,
            bounds=(0, None),
            units=pyunits.kg/pyunits.m**3,
            doc='Polymer solution density [kg/m3]')
        self.polymer_density.fix(1000)  # Low % in injection solution so assumed to be close to water

        self.polymer_ratio_in_soln = Var(initialize=0.03,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc='Polymer ratio in solution')
        self.polymer_ratio_in_soln.fix(0.003)

        self.polymer_capital_A = Var(initialize=30507,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Polymer capital A factor')
        self.polymer_capital_A.fix(30507)

        self.polymer_capital_B = Var(initialize=0.8274,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Polymer capital B factor')
        self.polymer_capital_B.fix(0.8274)

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

        self.polymer_feed_rate = Var(initialize=500,
            bounds=(0, None),
            units=pyunits.lb/pyunits.day,
            doc='Polymer feed rate [lb/day]')

        self.polymer_soln_flow= Var(initialize=500,
            bounds=(0, None),
            units=pyunits.gallon/pyunits.min,
            doc='Polymer solution feed rate [gal/min]')

        self.polymer_feed_rate_constr = Constraint(expr=
            self.polymer_feed_rate == pyunits.convert(
                self.flow_in * self.dose, 
                to_units=pyunits.lb/pyunits.day
            ))

        self.polymer_soln_flow_constr = Constraint(expr=
            self.polymer_soln_flow == pyunits.convert(pyunits.convert(
                self.polymer_feed_rate, to_units=pyunits.kg/pyunits.day
            ) / self.polymer_density / self.polymer_ratio_in_soln, 
            to_units=pyunits.gallon/pyunits.minute)
            )
        
        self.chem_dict = {'Anionic_Polymer': 0.5 * self.dose, 
                          'Cationic_Polymer': 0.5 * self.dose}


    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2007
        tpec_tic = 'TPEC'
        self.poly_setup()
        
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
            (self.polymer_capital_A * self.polymer_feed_rate ** self.polymer_capital_B) *
            self.tpec_tic * 1E-6,
            doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=
            ((0.746 * self.polymer_soln_flow * self.lift_height) / 
            (3960 * self.pump_eff * self.motor_eff)) / self.flow_in, 
            doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)