from pyomo.environ import Var, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.core.wt3_unit_sido import WT3UnitProcess

## REFERENCE
## CAPITAL:
# Poseidon
## ELECTRICITY:
# none

module_name = 'grit_chamber'

class UnitProcess(WT3UnitProcess):

    def grit_setup(self):
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.day)

        self.grit_cap_A_factor = Var(initialize=9.13003E-3,
            bounds=(0, None),
            doc='Grit chamber capital A factor')
        self.grit_cap_A_factor.fix(9.13003E-3)

        self.grit_cap_B_factor = Var(initialize=0.446445,
            bounds=(0, None),
            doc='Wire screen capital B factor')
        self.grit_cap_B_factor.fix(0.446445)

        self.grit_energy_A_factor = Var(initialize=4.135609,
            bounds=(0, None),
            doc='Grit chamber annual energy A factor')
        self.grit_energy_A_factor.fix(4.135609)

        self.grit_energy_B_factor = Var(initialize=1.007629,
            bounds=(0, None),
            doc='Wire screen annual energy B factor')
        self.grit_energy_B_factor.fix(1.007629)

        self.grit_other_A_factor = Var(initialize=0.900323E-3,
            bounds=(0, None),
            doc='Grit chamber annual energy A factor')
        self.grit_other_A_factor.fix(0.900323E-3)

        self.grit_other_B_factor = Var(initialize=0.443285,
            bounds=(0, None),
            doc='Wire screen annual energy B factor')
        self.grit_other_B_factor.fix(0.443285)

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2006
        tpec_tic = 'TIC'
        self.grit_setup()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
                (self.grit_cap_A_factor * self.flow_in ** self.grit_cap_B_factor) *
                self.tpec_tic,
                doc='Unadjusted fixed capital investment')
        self.costing.other_var_cost = \
            self.grit_other_A_factor * self.flow_in ** self.grit_other_B_factor
        self.electricity = Expression(expr=
                (self.grit_energy_A_factor * self.flow_in ** self.grit_energy_B_factor) /
                pyunits.convert(self.flow_in, to_units=pyunits.m**3/pyunits.year),
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)