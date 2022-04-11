from pyomo.environ import Var, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# Voutchkov Fig 4.7
## ELECTRICITY:
# none

module_name = 'micro_screen'

class UnitProcess(WT3UnitProcess):

    def micro_setup(self):
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.day)

        self.micro_screen_A_factor = Var(initialize=1.551E-5,
            bounds=(0, None),
            doc='Micro screen capital A factor')
        self.micro_screen_A_factor.fix(1.551E-5)

        self.micro_screen_B_factor = Var(initialize=0.959,
            bounds=(0, None),
            doc='Micro screen capital B factor')
        self.micro_screen_B_factor.fix(0.959)

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2018
        tpec_tic = 'TIC'
        self.micro_setup()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
                (self.micro_screen_A_factor * self.flow_in ** self.micro_screen_B_factor) *
                self.tpec_tic,
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=0,
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)