from pyomo.environ import Var, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# Voutchkov Fig 4.5
## ELECTRICITY:
# none

module_name = 'band_screen'

class UnitProcess(WT3UnitProcess):

    def band_setup(self):
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.day)

        self.band_screen_A_factor = Var(initialize=6.79E-6,
            bounds=(0, None),
            doc='Band screen capital A factor')
        self.band_screen_A_factor.fix(9.52E-6)

        self.band_screen_B_factor = Var(initialize=1.03,
            bounds=(0, None),
            doc='Band screen capital B factor')
        self.band_screen_B_factor.fix(1.03)

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2018
        tpec_tic = 'TIC'
        self.band_setup()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
                (self.band_screen_A_factor * self.flow_in ** self.band_screen_B_factor) *
                self.tpec_tic,
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=0,
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)