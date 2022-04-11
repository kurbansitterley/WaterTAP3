from pyomo.environ import Var, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# Voutchkov Fig 4.5
## ELECTRICITY:
# none

module_name = 'drum_screen'

class UnitProcess(WT3UnitProcess):

    def drum_setup(self):
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.day)

        self.drum_screen_A_factor = Var(initialize=1.59E-5,
            bounds=(0, None),
            doc='Drum screen capital A factor')
        self.drum_screen_A_factor.fix(1.59E-5)

        self.drum_screen_B_factor = Var(initialize=0.984,
            bounds=(0, None),
            doc='Drum screen capital B factor')
        self.drum_screen_B_factor.fix(0.984)

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2018
        tpec_tic = 'TIC'
        self.drum_setup()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
                (self.drum_screen_A_factor * self.flow_in ** self.drum_screen_B_factor) *
                self.tpec_tic,
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=0,
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)