from pyomo.environ import Constraint, Expression, Var, units as pyunits
from watertap3.utils import financials
from watertap3.core.wt3_unit_siso import WT3UnitProcessSISO

## REFERENCE
## CAPITAL:
# McGiveney & Kawamura
## ELECTRICITY:
# 

module_name = 'flocculator'

class UnitProcess(WT3UnitProcessSISO):

    def floc_setup(self):
        '''

        :return:
        '''
        time = self.flowsheet().config.time
        t = time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[t], 
            to_units=pyunits.m**3/pyunits.hr)
        
        self.residence_time = Var(
            initialize=10, 
            units=pyunits.minute, 
            bounds=(5,45), 
            doc='Flocculator residence time [min]')
        self.residence_time.fix(10)

        self.floc_motor_eff = Var(
            initialize=0.75,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Flocculator mixer motor efficiency')
        self.floc_motor_eff.fix(0.75)

        self.num_mixers = Var(
            initialize=3,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Number flocculator mixers')
        self.num_mixers.fix(3)

        self.g = Var(
            initialize=20,
            bounds=(0, None),
            units=pyunits.s**-1,
            doc='Velocity gradient')
        self.g.fix(20)

        self.viscosity = Var(
            initialize=1E-3,
            bounds=(0, None),
            units=pyunits.kilogram/(pyunits.second*pyunits.meter),
            doc='Water viscosity [kg/m*s]')
        self.viscosity.fix(1E-3)

        self.floc_basin_vol = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.Mgallons,
            doc='Flocculator basin volume [Mgal]')

        self.floc_capital_A = Var(
            initialize=1E6,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Flocculator capital A parameter')

        self.floc_capital_B = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc='Flocculator capital B parameter')

        self.floc_mixer_power = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.kW,
            doc='Flocculator mixer power [kW]')

        for k, v in self.unit_params.items():
            if k in ['residence_time', 'floc_motor_eff', 'num_mixers']:
                getattr(self, k).fix(v)
        
        if 'terminal_floc' in self.unit_params.keys():
            self.terminal_floc = self.unit_params['terminal_floc']
        else:
            self.terminal_floc = False

        if 'vel_gradient' in self.unit_params.keys():
            self.vel_gradient = self.unit_params['vel_gradient']
            if self.vel_gradient not in [20, 50, 80]:
                self.vel_gradient = 20
        else:
            self.vel_gradient = 20
        
        self.g.fix(self.vel_gradient)

        # if self.terminal_floc:
        #     self.removal_fraction[0, 'toc'].fix(0.40)

        self.floc_basin_vol_constr = Constraint(expr=self.floc_basin_vol ==
                pyunits.convert(self.flow_in, 
                to_units=pyunits.Mgallon/pyunits.minute) * self.residence_time)

        self.floc_mixer_power_constr = Constraint(expr=
            self.floc_mixer_power == pyunits.convert((self.g ** 2 * 
            pyunits.convert(self.floc_basin_vol,
                to_units=pyunits.m**3) * self.viscosity) * self.num_mixers, 
                to_units=pyunits.kilowatt) / self.floc_motor_eff)

        if self.vel_gradient == 20:
            self.floc_capital_A.fix(908675)
            self.floc_capital_B.fix(0.7191)
        if self.vel_gradient == 50:
            self.floc_capital_A.fix(1000977)
            self.floc_capital_B.fix(0.7579)
        if self.vel_gradient == 80:
            self.floc_capital_A.fix(1218085)
            self.floc_capital_B.fix(0.8266)

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2007
        tpec_tic = 'TPEC'
        self.floc_setup()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
                (self.floc_capital_A * self.floc_basin_vol ** self.floc_capital_B) *
                self.tpec_tic * 1E-6,
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.floc_mixer_power / self.flow_in,
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)

