from pyomo.environ import Block, Expression, Var, Constraint, NonNegativeReals, units as pyunits
from watertap3.utils import financials
from watertap3.wt_units.wt_unit import WT3UnitProcess

## REFERENCE
## CAPITAL:
# McGivney & Kawamura (2008) 
## ELECTRICITY:
# Biosolids Treatment Processes (2007)
# https://doi.org/10.1007/978-1-59259-996-7

## FACT SHEETS
# https://www.epa.gov/sites/default/files/2018-11/documents/recessed-plate-filter-press-factsheet.pdf
# https://www.epa.gov/sites/default/files/2018-11/documents/belt-filter-press-factsheet.pdf


module_name = 'filter_press'
basis_year = 2007
tpec_or_tic = 'TPEC'



class UnitProcess(WT3UnitProcess):

    def fp_setup(self):

        time = self.flowsheet().config.time
        t = time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[t], 
            to_units=pyunits.m ** 3 / pyunits.hr)
        self.chem_dict = {}

        self.hours_per_day_operation = Var(initialize=24,
            units=pyunits.hours/pyunits.day,
            bounds=(0, 24.1),
            doc='Hours per day of operation [hr/day]')
        self.hours_per_day_operation.fix(24)

        self.cycle_time = Var(initialize=2,
            units=pyunits.hours,
            bounds=(0, None),
            doc='Filter press cycle time [hr]')
        self.cycle_time.fix(3)

        self.cycles_per_day = Var(initialize=8,
            units=pyunits.day**-1,
            bounds=(0, None),
            doc='Number of cycles per day [1/day]')

        self.press_capacity = Var(initialize=1000,
            units=pyunits.ft**3,
            bounds=(0, None),
            doc='Filter press capacity [ft3]')

        self.fp_capital_A = Var(initialize=8060.2,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc='Filter press capital A factor')

        self.fp_capital_B = Var(initialize=0.6015,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc='Filter press capital B factor')

        self.fp_electricity_A = Var(initialize=16.285,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc='Filter press electricity A factor')

        self.fp_electricity_B = Var(initialize=1.2434,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc='Filter press electricity B factor')

        for k, v in self.unit_params.items():
            if k in ['hours_per_day_operation', 'cycle_time']:
                getattr(self, k).fix(v)
        
        if 'type' in self.unit_params.keys():
            self.filter_press_type = self.unit_params['type']
            if self.filter_press_type not in ['belt', 'pressure']:
                self.filter_press_type = 'belt'
        else:
            self.filter_press_type = 'belt'

        self.cycles_per_day_constr = Constraint(expr=self.cycles_per_day == 
            self.hours_per_day_operation / self.cycle_time)

        self.press_capacity_constr = Constraint(expr=self.press_capacity == 
            pyunits.convert(self.flow_in, to_units=pyunits.ft**3/pyunits.day) /
            self.cycles_per_day)

        if self.filter_press_type == 'belt':
            self.fp_capital_A.fix(8060.2)
            self.fp_capital_B.fix(0.6015)
            self.fp_electricity_A.fix(16.285)
            self.fp_electricity_B.fix(1.2434)

        if self.filter_press_type == 'pressure':
            self.fp_capital_A.fix(102794)
            self.fp_capital_B.fix(0.4216)
            self.fp_electricity_A.fix(16.612)
            self.fp_electricity_B.fix(1.2195)

    def get_costing(self, unit_params=None, year=None):
        '''
        Initialize the unit in WaterTAP3.
        '''
        financials.create_costing_block(self, basis_year, tpec_or_tic)
        self.fp_setup()
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=
            (self.fp_capital_A * pyunits.convert(self.flow_in, 
            to_units=pyunits.gallon/pyunits.hour) ** self.fp_capital_B) * 1E-6,
            doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=
            (self.fp_electricity_A * self.press_capacity ** self.fp_electricity_B) / 
            pyunits.convert(self.flow_in, to_units=pyunits.m**3/pyunits.year),
            doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing)
