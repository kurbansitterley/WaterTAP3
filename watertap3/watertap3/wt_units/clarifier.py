from pyomo.environ import Var, Constraint, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# Cost Estimating Manual for Water Treatment Facilities (McGivney/Kawamura) (2008)
# DOI:10.1002/9780470260036
# Water and Wastewater Engineering: Design Principles and Practice (Mackenzie L. Davis) (2010)
## ELECTRICITY:


module_name = 'clarifier'
basis_year = 2007
tpec_or_tic = 'TIC'


class UnitProcess(WT3UnitProcess):

    def sed_setup(self):
        time = self.flowsheet().config.time.first()
        self.chem_dict = {}
        self.residence_time = Var(
            initialize=2, 
            units=pyunits.hr, 
            bounds=(0, None), 
            doc='Residence time [hr]')
        self.residence_time.fix(2)

        self.basin_height = Var(
            initialize=10, 
            units=pyunits.ft, 
            bounds=(0, None), 
            doc='Basin height [ft]')
        self.basin_height.fix(10)
    
        self.clarifier_capital_A = Var(
            initialize=3470.6,
            units=pyunits.dimensionless, 
            bounds=(0, None), 
            doc='Clarifier capital A factor')
        self.clarifier_capital_A.fix(3470.6)

        self.clarifier_capital_B = Var(
            initialize=0.6173, 
            units=pyunits.dimensionless,  
            bounds=(0, None), 
            doc='Clarifier capital B factor')
        self.clarifier_capital_B.fix(0.6173)

        self.basin_surface_area = Var(
            initialize=10000, 
            units=pyunits.ft**2, 
            bounds=(0, None), 
            doc='Basin surface area [ft2]')

        for k, v in self.unit_params.items():
            if k in ['residence_time', 'basin_height']:
                getattr(self, k).fix(v)
        if 'water_recovery' in self.unit_params.keys():
            self.water_recovery.fix(self.unit_params['water_recovery'])

        self.basin_surface_area_constr = Constraint(expr=
            self.basin_surface_area == pyunits.convert(
                (self.flow_vol_in[time] * self.residence_time) / self.basin_height,
                to_units=pyunits.ft**2
            ))

    def get_costing(self, unit_params=None, year=None):
        '''
        Initialize the unit in WaterTAP3.
        '''
        financials.create_costing_block(self, basis_year, tpec_or_tic)
        self.sed_setup()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
                (self.clarifier_capital_A * self.basin_surface_area ** self.clarifier_capital_B) *
                self.tpec_tic * 1E-6,
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=0,
                                      doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing)