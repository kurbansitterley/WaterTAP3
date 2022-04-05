from pyomo.environ import Constraint, Expression, Var, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# McGiveney & Kawamura
## ELECTRICITY:
# 

module_name = 'rapid_mix'

class UnitProcess(WT3UnitProcess):

    def rapid_mix_setup(self):
        '''
        :return:
        '''
        time = self.flowsheet().config.time
        t = time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[t], 
            to_units=pyunits.m**3/pyunits.hr)

        self.residence_time = Var(
            initialize=5, 
            units=pyunits.second, 
            bounds=(5, 60), 
            doc='Rapid mix residence time [s]')
        self.residence_time.fix(5)

        self.mixer_motor_eff = Var(
            initialize=0.75,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc='Rapid mixer motor efficiency')
        self.mixer_motor_eff.fix(0.75)

        self.g = Var(
            initialize=300,
            bounds=(0, None),
            units=pyunits.s**-1,
            doc='Velocity gradient')
        self.g.fix(300)

        self.viscosity = Var(
            initialize=1E-3,
            bounds=(0, None),
            units=pyunits.kilogram/(pyunits.second*pyunits.meter),
            doc='Water viscosity [kg/m*s]')
        self.viscosity.fix(1E-3)

        self.rapid_mix_basin_vol = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.gallons,
            doc='Rapid mixer basin volume [gal]')

        self.rm_capital_A = Var(
            initialize=1E6,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Rapid mixer capital A parameter')

        self.rm_capital_B = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Rapid mixer capital B parameter')

        self.rapid_mixer_power = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.kW,
            doc='Rapid mixer power [kW]')

        for k, v in self.unit_params.items():
            if k in ['residence_time']:
                getattr(self, k).fix(v)

        if 'vel_gradient' in self.unit_params.keys():
            self.vel_gradient = self.unit_params['vel_gradient']
            if self.vel_gradient not in [300, 600, 900]:
                self.vel_gradient = 300
        else:
            self.vel_gradient = 300
        
        self.g.fix(self.vel_gradient)

        self.rapid_mix_basin_vol_constr = Constraint(expr=self.rapid_mix_basin_vol ==
                pyunits.convert(self.flow_in * self.residence_time, 
                to_units=pyunits.gallons))

        self.rapid_mixer_power_constr = Constraint(expr=
            self.rapid_mixer_power == pyunits.convert((self.g ** 2 * 
            pyunits.convert(self.rapid_mix_basin_vol, to_units=pyunits.m**3) 
            * self.viscosity), to_units=pyunits.kilowatt) / 
            self.mixer_motor_eff)

        if self.vel_gradient == 300:
            self.rm_capital_A.fix(1268.1)
            self.rm_capital_B.fix(0.4431)
        if self.vel_gradient == 600:
            self.rm_capital_A.fix(1131.2)
            self.rm_capital_B.fix(0.4731)
        if self.vel_gradient == 900:
            self.rm_capital_A.fix(618.01)
            self.rm_capital_B.fix(0.5696)

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2007
        tpec_tic = 'TPEC'
        self.rapid_mix_setup()
        
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
                (self.rm_capital_A * self.rapid_mix_basin_vol ** self.rm_capital_B) *
                self.tpec_tic * 1E-6,
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.rapid_mixer_power / self.flow_in,
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)
