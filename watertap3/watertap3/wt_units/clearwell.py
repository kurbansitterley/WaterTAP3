from pyomo.environ import Var, Constraint, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# citation here
## ELECTRICITY:
# citation here

module_name = 'unit'
basis_year = 2020
tpec_or_tic = 'TPEC'


class UnitProcess(WT3UnitProcess):

    def clearwell_setup(self):
        time = self.flowsheet().config.time.first()
        self.chem_dict = {}
        self.flow_in = pyunits.convert(self.flow_vol_in[time], 
            to_units=pyunits.Mgallons/pyunits.hour)

        self.storage_time = Var(
            initialize=3,
            bounds=(0, None),
            units=pyunits.hour,
            doc='Storage time for clearwell [hr]')
        self.storage_time.fix(3)

        self.storage_vol = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.Mgallons,
            doc='Storage volume of clearwell [MG]')

        self.clearwell_capital_A = Var(
            initialize=1016775,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Clearwell capital A factor')
        self.clearwell_capital_A.fix(1016775)

        self.clearwell_capital_B = Var(
            initialize=0.6836,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Clearwell capital B factor')
        self.clearwell_capital_B.fix(0.6836)

        if 'storage_time' in self.unit_params.keys():
            self.storage_time.fix(self.unit_params['storage_time'])

        self.storage_vol_constr = Constraint(expr=
            self.storage_vol == self.flow_in * self.storage_time)


    def get_costing(self, unit_params=None, year=None):
        '''
        Initialize the unit in WaterTAP3.
        '''
        self.clearwell_setup()
        financials.create_costing_block(self, basis_year, tpec_or_tic)
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
            (self.clearwell_capital_A * self.storage_vol ** self.clearwell_capital_B) * 1E-6,
            doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=0,
                                      doc='Electricity intensity [kwh/m3]')
        financials.get_complete_costing(self.costing)