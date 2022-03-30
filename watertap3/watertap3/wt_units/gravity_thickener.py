from pyomo.environ import Var, Constraint, Expression, sqrt, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# citation here
## ELECTRICITY:
# citation here

module_name = 'gravity_thickener'
basis_year = 2007
tpec_or_tic = 'TPEC'


class UnitProcess(WT3UnitProcess):

    def grav_thick_setup(self):
        t = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[t], 
            to_units=pyunits.m**3/pyunits.hr)

        self.chem_dict = {}

        self.pct_solids = Var(initialize=0.06,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Percent influent solids [lb/lb]')
        self.pct_solids.fix(0.06)

        self.density = Var(initialize=1000,
            bounds=(0, None),
            units=pyunits.kg/pyunits.m**3,
            doc='Influent density [kg/m3]')

        self.mass_loading = Var(initialize=1,
            bounds=(0.1, 1.5),
            units=pyunits.lb/pyunits.ft**2/pyunits.day,
            doc='Mass loading [lb/ft2/d]')
        self.mass_loading.fix(0.7)
        
        self.grav_thick_capital_A = Var(initialize=2798.7,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Gravity thickener capital A')
        self.grav_thick_capital_A.fix(2798.7)

        self.grav_thick_capital_B = Var(initialize=1.305,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Gravity thickener capital B')
        self.grav_thick_capital_B.fix(1.305)

        self.grav_thick_energy_A = Var(initialize=15.077,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Gravity thickener energy A')
        self.grav_thick_energy_A.fix(15.077)

        self.grav_thick_energy_B = Var(initialize=0.8287,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Gravity thickener energy B')
        self.grav_thick_energy_B.fix(0.8287)

        self.area_required = Var(initialize=1000,
            bounds=(0, None),
            units=pyunits.ft**2,
            doc='Area required per gravity thickener [ft2]')

        self.num_thickeners = Var(initialize=3,
            bounds=(1, 10),
            units=pyunits.dimensionless,
            doc='Number gravity thickeners')
        self.num_thickeners.fix(4)

        self.diameter = Var(initialize=1,
            bounds=(0, None),
            units=pyunits.ft,
            doc='Thickener diameter [ft]')

        self.density_constr = Constraint(expr=
            self.density == 0.6312 * 
            sum([self.conc_mass_in[t, c] for c in self.config.property_package.component_list]) + 997.86)
        
        self.area_constr = Constraint(expr=(self.area_required * self.num_thickeners) == 
            (pyunits.convert(self.flow_in, to_units=pyunits.gallon/pyunits.day) * 
            self.pct_solids * 
            pyunits.convert(self.density, to_units=pyunits.lb/pyunits.gallon)) / 
            self.mass_loading)

        self.diameter_constr = Constraint(expr=self.diameter == 
            sqrt((self.area_required * 4) / 3.14159))

    def get_costing(self, unit_params=None, year=None):
        '''
        Initialize the unit in WaterTAP3.
        '''
        self.grav_thick_setup()
        financials.create_costing_block(self, basis_year, tpec_or_tic)
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
            (self.grav_thick_capital_A * self.diameter ** self.grav_thick_capital_B) * 
            self.num_thickeners * 1E-6,
            doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=
            (self.grav_thick_energy_A * self.area_required ** self.grav_thick_energy_B * self.num_thickeners) /
            pyunits.convert(self.flow_in, to_units=pyunits.m**3/pyunits.year),
            doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing)