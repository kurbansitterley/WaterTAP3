from pyomo.environ import Var, Constraint, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.core.wt3_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# Based on costs for FILTER MEDIA DUAL MEDIA - FIGURE 5.5.27
# Cost Estimating Manual for Water Treatment Facilities (McGivney/Kawamura) (2008)
# DOI:10.1002/9780470260036
## ELECTRICITY:
# An Analysis of Energy Consumption and the Use of Renewables for a Small Drinking Water Treatment Plant
# Water 2020, 12, 28; doi:10.3390/w12010028

module_name = 'media_filtration'

class UnitProcess(WT3UnitProcess):

    def media_filter_setup(self):

        time = self.flowsheet().config.time.first()
        
        self.flow_in = pyunits.convert(self.flow_vol_in[time], 
            to_units=pyunits.m**3/pyunits.hr)

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

        self.media_capital_A = Var(initialize=1E3,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc='Media capital A factor')

        self.media_capital_B = Var(initialize=0.7,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc='Media capital B factor')

        self.structure_capital_A = Var(initialize=10439,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc='Filter structure capital A factor')
        self.structure_capital_A.fix(10439)

        self.structure_capital_B = Var(initialize=0.6922,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc='Filter structure capital B factor')
        self.structure_capital_B.fix(0.6922)

        self.backwash_capital_A = Var(initialize=10439,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc='Backwash capital A factor')
        self.backwash_capital_A.fix(10056)

        self.backwash_capital_B = Var(initialize=0.5222,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc='Backwash capital B factor')
        self.backwash_capital_B.fix(0.5222)

        self.filter_surface_area = Var(initialize=100,
            units=pyunits.ft**2,
            bounds=(0, None),
            doc='Surface area [ft2]')

        self.media_cost = Var(initialize=1E3,
            bounds=(0, None),
            doc='Media cost')

        self.filter_backwash_cost = Var(initialize=1E3,
            bounds=(0, None),
            doc='Backwash cost')

        self.filter_structure_cost = Var(initialize=1E3,
            bounds=(0, None),
            doc='Structure cost')

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

        if 'media_type' in self.unit_params.keys():
            self.media_type = self.unit_params['media_type']
            if self.media_type not in ['dual', 'tri', 'sand']:
                self.media_type == 'dual'
        else:
            self.media_type = 'dual'

        if self.media_type == 'dual':
            self.media_capital_A.fix(98.0182)
            self.media_capital_B.fix(0.9056)
        if self.media_type == 'tri':
            self.media_capital_A.fix(207.97)
            self.media_capital_B.fix(0.8729)
        if self.media_type == 'sand':
            self.media_capital_A.fix(296.15)
            self.media_capital_B.fix(0.714)

        self.surface_area_constr = Constraint(expr=
            self.filter_surface_area == pyunits.convert(self.flow_in / 
            (self.loading_rate * self.number_units),
            to_units=pyunits.ft**2))

        self.media_capital_constr = Constraint(expr=
            self.media_cost == (self.media_capital_A * self.filter_surface_area ** self.media_capital_B) * 1E-6)

        self.backwash_cost_constr = Constraint(expr=
            self.filter_backwash_cost == 
            ((self.backwash_capital_A * self.filter_surface_area ** self.backwash_capital_B) * 
            (1 - self.water_recovery[time])) * 1E-6)

        self.structure_cost_constr = Constraint(expr=
            self.filter_structure_cost == (self.structure_capital_A * \
                self.filter_surface_area ** self.structure_capital_B) * 1E-6)

        self.media_filter_fixed_cap_constr = Constraint(expr=
            self.media_filter_fixed_cap == #self.tpec_tic * \
            ((self.media_cost + self.filter_structure_cost + self.filter_backwash_cost) 
            * self.number_units))


    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2007
        tpec_tic = 'TPEC'
        self.media_filter_setup()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.media_filter_fixed_cap,
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.media_filter_electricity,
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)