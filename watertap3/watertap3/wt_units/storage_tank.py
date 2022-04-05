from pyomo.environ import Var, Constraint, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# "Cone roof tank" costs from:
# DOE/NETL-2002/1169 - Process Equipment Cost Estimation Final Report
# Loh, H. P., Lyons, Jennifer, and White, Charles W. Process Equipment Cost Estimation, Final Report.
# United States: N. p., 2002. Web. doi:10.2172/797810.
# Regression of cost vs. capacity
# Capacity calculated based on storage time (user input)

module_name = 'storage_tank'
basis_year = 1998
tpec_or_tic = 'TPEC'

class UnitProcess(WT3UnitProcess):

    def tank_setup(self):
        time = self.flowsheet().config.time.first()
        self.chem_dict = {}
        self.flow_in = pyunits.convert(self.flow_vol_in[time], 
            to_units=pyunits.m**3/pyunits.hr)

        self.tank_capital_A = Var(initialize=0.00344,
            bounds=(0, None),
            doc='Storage tank capital A factor')
        self.tank_capital_A.fix(0.00344)

        self.tank_capital_B = Var(initialize=3,
            bounds=(0, None),
            doc='Storage tank capital B factor')
        self.tank_capital_B.fix(0.72093)

        self.storage_time = Var(initialize=6,
            bounds=(0, None),
            units=pyunits.hour,
            doc='Storage duration [hour]')
        self.storage_time.fix(6)

        self.surge_capacity = Var(initialize=0.2,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Storage tank surge capacity [%]')
        self.surge_capacity.fix(0.2)

        self.storage_vol = Var(initialize=1000,
            bounds=(0, None),
            units=pyunits.m**3,
            doc='Storage tank volume [m3]')

        for k, v in self.unit_params.items():
            if k in ['storage_time', 'surge_capacity']:
                getattr(self, k).fix(v)
        
        self.storage_vol_constr = Constraint(expr= 
            self.storage_vol == self.flow_in * self.storage_time *
                (1 + self.surge_capacity))

    def get_costing(self, unit_params=None, year=None):
        '''
        Initialize the unit in WaterTAP3.
        '''
        self.tank_setup()
        financials.create_costing_block(self, basis_year, tpec_or_tic)
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
            self.tank_capital_A * self.storage_vol ** self.tank_capital_B ,
            doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=0,
            doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing)