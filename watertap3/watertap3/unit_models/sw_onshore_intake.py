from pyomo.environ import Expression, units as pyunits
from watertap3.utils import financials
from watertap3.core.wt3_unit_pt import WT3UnitProcessPT

## REFERENCE: Voutchkov (2018) figures 4.2 and 4.4

module_name = 'sw_onshore_intake'

class UnitProcess(WT3UnitProcessPT):

    def fixed_cap(self):
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3 / pyunits.day)
        self.base_fixed_cap_cost = 0.000215
        self.cap_scaling_exp = 0.888803
        sw_intake_cap = self.base_fixed_cap_cost * self.flow_in ** self.cap_scaling_exp
        return sw_intake_cap

    def elect(self):
        self.lift_height = 100 * pyunits.ft
        self.pump_eff = 0.9 * pyunits.dimensionless
        self.motor_eff = 0.9 * pyunits.dimensionless
        flow_in_gpm = pyunits.convert(self.flow_in,
            to_units=pyunits.gallons/pyunits.minute)
        flow_in_m3hr = pyunits.convert(self.flow_in,
            to_units=pyunits.m**3 / pyunits.hour)
        electricity = (0.746 * flow_in_gpm * self.lift_height / \
            (3960 * self.pump_eff * self.motor_eff)) / flow_in_m3hr
        return electricity

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2018
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.elect(),
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year)


