from pyomo.environ import Var, Constraint, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# McGivney/Kawamura (2008) Fig. 5.5.16
## ELECTRICITY:
# citation here

module_name = 'pac_addition'

class UnitProcess(WT3UnitProcess):

    def pac_setup(self):

        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time], 
            to_units=pyunits.m**3/pyunits.hr)

        self.pac_dose = Var(initialize=1,
            bounds=(0, None),
            units=pyunits.mg/pyunits.liter,
            doc='PAC dose [mg/L]')
        self.pac_dose.fix(10)

        if 'dose' in self.unit_params.keys():
            self.pac_dose.fix(self.unit_params['dose'])

        self.pac_capital_A = Var(initialize=45190,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='PAC capital A factor')
        self.pac_capital_A.fix(45190)

        self.pac_capital_B = Var(initialize=0.3233,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='PAC capital B factor')
        self.pac_capital_B.fix(0.3233)

        self.pac_feed_rate = Var(initialize=500,
            bounds=(0, None),
            units=pyunits.lb/pyunits.hr,
            doc='PAC feed rate [lb/hr]')

        self.pac_feed_rate_constr = Constraint(expr=
            self.pac_feed_rate == pyunits.convert(
                self.flow_in * self.pac_dose, 
                to_units=pyunits.lb/pyunits.hr
            ))
        
        self.chem_dict = {'Powder_Activated_Carbon': 
            pyunits.convert(self.pac_dose, 
            to_units=pyunits.kg/pyunits.m**3)}


    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2007
        tpec_tic = 'TPEC'
        self.pac_setup()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
            (self.pac_capital_A * self.pac_feed_rate ** self.pac_capital_B) *
            self.tpec_tic * 1E-6,
            doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=0,
            doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)