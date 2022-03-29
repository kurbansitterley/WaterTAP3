from pyomo.environ import Var, Constraint, Expression, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# Cost Estimating Manual for Water Treatment Facilities (McGivney/Kawamura) (2008)
# DOI:10.1002/9780470260036
# Water and Wastewater Engineering: Design Principles and Practice (Mackenzie L. Davis) (2010)
## ELECTRICITY:


module_name = 'sedimentation'
basis_year = 2007
tpec_or_tic = 'TPEC'


class UnitProcess(WT3UnitProcess):

    def sed_setup(self):
        time = self.flowsheet().config.time.first()
        self.chem_dict = {}
        self.settling_velocity = Var(
            initialize=0.005, 
            units=pyunits.m/pyunits.second, 
            bounds=(0, None), 
            doc='Settling velocity [m/s]')
        self.settling_velocity.fix(0.005)
    
        self.sed_basin_capital_A = Var(
            initialize=13572, 
            units=pyunits.dimensionless, 
            bounds=(0, None), 
            doc='Sedimentation basin capital A factor')
        self.sed_basin_capital_A.fix(13572)

        self.sed_basin_capital_B = Var(
            initialize=0.3182, 
            units=pyunits.dimensionless,  
            bounds=(0, None), 
            doc='Sedimentation basin capital B factor')
        self.sed_basin_capital_B.fix(0.3182)

        self.basin_surface_area = Var(
            initialize=10000, 
            units=pyunits.ft**2, 
            bounds=(0, None), 
            doc='Basin surface area [ft2]')

        if 'settling_velocity' in self.unit_params.keys():
            self.settling_velocity.fix(self.unit_params['settling_velocity'])
        if 'water_recovery' in self.unit_params.keys():
            self.water_recovery.fix(self.unit_params['water_recovery'])
        else:
            self.water_recovery.fix(0.99)

        self.basin_surface_area_constr = Constraint(expr=
            self.basin_surface_area == pyunits.convert(
                self.flow_vol_in[time] / self.settling_velocity,
                to_units=pyunits.ft**2
            ))

    def get_costing(self, unit_params=None, year=None):
        '''
        Initialize the unit in WaterTAP3.
        '''
        financials.create_costing_block(self, basis_year, tpec_or_tic)
        self.sed_setup()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
                (self.sed_basin_capital_A * self.basin_surface_area ** self.sed_basin_capital_B) *
                self.tpec_tic * 1E-6,
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=0,
                                      doc='Electricity intensity [kwh/m3]')
        financials.get_complete_costing(self.costing)