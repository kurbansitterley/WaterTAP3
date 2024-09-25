from pyomo.environ import Expression, units as pyunits
from watertap3.utils import financials
from watertap3.core.wt3_unit_pt import WT3UnitProcessPT

## REFERENCE
## CAPITAL:
# Based on TABLE 7.2 in:
# Chemical Engineering Design, 2nd Edition. Principles, Practice and Economics of Plant and Process Design (2012)
# https://www.elsevier.com/books/chemical-engineering-design/towler/978-0-08-096659-5
# eBook ISBN: 9780080966601
## ELECTRICITY:

module_name = 'static_mixer'

class UnitProcess(WT3UnitProcessPT):

    def fixed_cap(self):
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.hr)
        self.number_of_units = 2
        self.a = 1065.7
        self.b = 0.336
        source_cost = self.a * self.flow_in ** self.b
        mix_cap = (source_cost * self.tpec_tic * self.number_of_units) * 1E-6
        return mix_cap

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2010
        tpec_tic = 'TPEC'
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=0,
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)