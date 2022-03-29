from pyomo.environ import Var, Constraint, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# Based on costs for FILTER MEDIA DUAL MEDIA - FIGURE 5.5.27
# Cost Estimating Manual for Water Treatment Facilities (McGivney/Kawamura) (2008)
# DOI:10.1002/9780470260036
## ELECTRICITY:
# An Analysis of Energy Consumption and the Use of Renewables for a Small Drinking Water Treatment Plant
# Water 2020, 12, 28; doi:10.3390/w12010028

module_name = 'media_filtration'
basis_year = 2007
tpec_or_tic = 'TPEC'


class UnitProcess(WT3UnitProcess):

    def media_filter_setup(self):

        time = self.flowsheet().config.time.first()
        self.chem_dict = {}
        self.flow_in = pyunits.convert(self.flow_vol_in[time], 
            to_units=pyunits.m ** 3 / pyunits.hr)

        self.loading_rate = Var(initialize=10,
            units=pyunits.m/pyunits.hr,
            bounds=(0, None),
            doc='Loading rate [m/hr]')
        self.loading_rate.fix(10)

        self.number_units = Var(initialize=6,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc='Number of units')
        self.number_units.fix(6)

        self.filter_surface_area = Var(initialize=100,
            units=pyunits.ft**2,
            bounds=(0, None),
            doc='Surface area [ft2]')

        self.dual_filter_cost = Var(initialize=1E3,
            bounds=(0, None),
            doc='Dual media filter cost')

        self.filter_backwash_cost = Var(initialize=1E3,
            bounds=(0, None),
            doc='Backwash cost')

        self.media_filter_fixed_cap = Var(initialize=1000,
            bounds=(0, None),
            doc='Media filtration capital [$M]')

        self.media_filter_electricity = Var(initialize=1.5E-4,
            units=pyunits.kWh/pyunits.m**3,
            bounds=(0, None),
            doc='Media filtration electricity intensity [kWh/m3]')
        self.media_filter_electricity.fix(1.5E-4)

        for k, v in self.unit_params.items():
            if k in ['loading_rate', 'number_units']:
                getattr(self, k).fix(v)

        self.surface_area_constr = Constraint(expr=
            self.filter_surface_area == pyunits.convert(self.flow_in / (self.loading_rate * self.number_units),
            to_units=pyunits.ft**2))

        self.dual_filter_cost_constr = Constraint(expr=
            self.dual_filter_cost == (38.319 * self.filter_surface_area + 21377))

        self.backwash_cost_constr = Constraint(expr=
            self.filter_backwash_cost == 292.44 * self.filter_surface_area + 92497)

        self.media_filter_fixed_cap_constr = Constraint(expr=
            self.media_filter_fixed_cap == self.tpec_tic * \
            ((self.dual_filter_cost + self.filter_backwash_cost) * self.number_units) * 1E-6)


    def get_costing(self, unit_params=None, year=None):
        '''
        Initialize the unit in WaterTAP3.
        '''
        financials.create_costing_block(self, basis_year, tpec_or_tic)
        self.media_filter_setup()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.media_filter_fixed_cap,
                                                           doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.media_filter_electricity,
                                      doc='Electricity intensity [kwh/m3]')
        financials.get_complete_costing(self.costing)