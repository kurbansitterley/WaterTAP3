from pyomo.environ import Var, Constraint, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.core.wt3_unit_sido import WT3UnitProcess

## REFERENCE
## CAPITAL:
# citation here
## ELECTRICITY:
# citation here

module_name = 'unit'

class UnitProcess(WT3UnitProcess):

    def fixed_cap(self):
        '''
        Docstrings go here.

        :return:
        '''
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.hr)
        unit_cap = 0
        return unit_cap

    def elect(self):
        '''
        Docstrings go here.

        :return:
        '''
        electricity = 0
        return electricity

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2020
        tpec_tic = 'TPEC'
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.elect(),
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)