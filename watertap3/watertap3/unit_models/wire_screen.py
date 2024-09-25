from pyomo.environ import Var, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.core.wt3_unit_sido import WT3UnitProcess

## REFERENCE
## CAPITAL:
# Voutchkov Fig 4.6
## ELECTRICITY:
# none

module_name = 'wire_screen'

class UnitProcess(WT3UnitProcess):

    def wire_setup(self):
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.day)

        self.wire_screen_A_factor = Var(initialize=9.52E-6,
            bounds=(0, None),
            doc='Wire screen capital A factor')
        self.wire_screen_A_factor.fix(9.52E-6)

        self.wire_screen_B_factor = Var(initialize=0.97457,
            bounds=(0, None),
            doc='Wire screen capital B factor')
        self.wire_screen_B_factor.fix(0.97457)

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2018
        tpec_tic = 'TIC'
        self.wire_setup()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
                (self.wire_screen_A_factor * self.flow_in ** self.wire_screen_B_factor) *
                self.tpec_tic,
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=0,
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)