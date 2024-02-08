from multiprocessing.sharedctypes import Value
from pyomo.environ import Expression, Var, units as pyunits
from watertap3.utils import financials
from watertap3.core.wt3_unit_sido import WT3UnitProcess

## REFERENCE: McGiveney & Kawamura

module_name = 'coag_and_floc'

class UnitProcess(WT3UnitProcess):

    def fixed_cap(self):
        '''

        :param unit_params: Input parameters from input sheet
        :type unit_params: dict

        :return:
        '''
        time = self.flowsheet().config.time
        t = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[t], 
            to_units=pyunits.m**3/pyunits.hr)
        flow_in_gpm = pyunits.convert(self.flow_in, 
            to_units=pyunits.gallons/pyunits.minute)
        self.alum_dose = Var(
            initialize=0.01, 
            units=pyunits.kg/pyunits.m**3, 
            doc='Alum dose [kg/m3]')
        self.alum_dose.fix(0.01)
        self.polymer_dose = Var( 
            initialize=0.0001, 
            units=pyunits.kg/pyunits.m**3, 
            doc='Polymer dose [kg/m3]')
        self.polymer_dose.fix(0.0001)
        if 'alum_dose' in self.unit_params.keys():
            self.alum_dose.fix(self.unit_params['alum_dose'] * 1E-3)
            # alum_dose_kgm3 = pyunits.convert(self.unit_params['alum_dose'] * (pyunits.mg/pyunits.liter), 
            #     to_units=(pyunits.kg/pyunits.m**3))
        # self.alum_dose.fix(alum_dose_kgm3)
        if 'polymer_dose' in self.unit_params.keys():
            self.polymer_dose.fix(self.unit_params['polymer_dose'] * 1E-3)
        # self.polymer_dose = pyunits.convert(self.unit_params['polymer_dose'] * (pyunits.mg/pyunits.liter), 
        #     to_units=(pyunits.kg/pyunits.m**3))
        # self.an_polymer = self.polymer_dose / 2
        # self.cat_polymer = self.polymer_dose / 2
        self.chem_dict = {'Aluminum_Al2_SO4_3': self.alum_dose, 
                'Anionic_Polymer': self.polymer_dose * 0.5, 
                'Cationic_Polymer': self.polymer_dose * 0.5}

        self.rapid_mixers = 1 * pyunits.dimensionless
        self.floc_mixers = 3 * pyunits.dimensionless
        self.rapid_mix_processes = 1 * pyunits.dimensionless
        self.floc_processes = 2 * pyunits.dimensionless
        self.coag_processes = 1 * pyunits.dimensionless
        self.floc_injection_processes = 1 * pyunits.dimensionless
        # self.rapid_mix_retention_time = 5.5 * pyunits.seconds
        # self.floc_retention_time = 12 * pyunits.minutes

        self.rapid_mix_retention_time = Var(time, 
            initialize=5.5, 
            units=pyunits.seconds, 
            doc='Rapid Mix Retention Time [s]')
        self.rapid_mix_retention_time.fix(5.5)
        self.floc_retention_time = Var(time, 
            initialize=12, 
            units=pyunits.minutes, 
            doc='Floc Retention Time [min]')
        self.floc_retention_time.fix(12)

        for k, v in self.unit_params.items():
            if k in ['rapid_mix_retention_time', 'floc_mix_retention_time']:
                getattr(self, k).fix(v)
        
        self.rapid_mix_basin_volume = pyunits.convert(self.rapid_mix_retention_time[t],
            to_units=pyunits.minutes) * flow_in_gpm
        self.floc_basin_volume = pyunits.convert(self.floc_retention_time[t] * flow_in_gpm,
            to_units=pyunits.Mgallons)
        self.alum_flow = self.alum_dose * self.flow_in
        self.alum_flow = pyunits.convert(self.alum_flow,
            to_units=(pyunits.lb/pyunits.hour))
        self.poly_flow = self.polymer_dose * self.flow_in
        self.poly_flow = pyunits.convert(self.poly_flow,
            to_units=(pyunits.lb/pyunits.day))
        self.rapid_mix_cap = (7.0814 * self.rapid_mix_basin_volume + 33269) * self.rapid_mix_processes # Rapid Mix, G = 900
        self.floc_cap = (952902 * self.floc_basin_volume + 177335) * self.floc_processes # Flocculator, G = 80
        self.coag_inj_cap = (212.32 * self.alum_flow + 73225) * self.coag_processes
        self.floc_inj_cap = (13662 * self.poly_flow + 20861) * self.floc_injection_processes
        self.coag_floc_cap = (self.rapid_mix_cap + self.floc_cap + \
            self.coag_inj_cap + self.floc_inj_cap) * 1E-6 * self.tpec_tic
        return self.coag_floc_cap

    def elect(self):
        '''
        Electricity intensity for coagulation/flocculation [kWh/m3]

        :return:
        '''
        flow_in_gpm = pyunits.convert(self.flow_in,
            to_units=pyunits.gallons/pyunits.minute)
        self.rapid_mix_basin_volume = pyunits.convert(self.rapid_mix_retention_time[0],
            to_units=pyunits.minutes) * flow_in_gpm
        self.rapid_mix_basin_volume = pyunits.convert(self.rapid_mix_basin_volume,
            to_units=pyunits.m**3)
        self.rapid_mix_power_consumption = (900 ** 2 * 0.001 * \
            self.rapid_mix_basin_volume) * (pyunits.watts/pyunits.m**3)
        self.rapid_mix_power = pyunits.convert(self.rapid_mix_power_consumption * self.rapid_mixers,
            to_units=pyunits.kilowatts)
        self.floc_basin_volume = self.floc_retention_time[0] * flow_in_gpm
        self.floc_basin_volume = pyunits.convert(self.floc_basin_volume,
            to_units=(pyunits.m**3))
        self.floc_power_consumption = (80 ** 2 * 0.001 * \
            self.floc_basin_volume) * (pyunits.watts/pyunits.m**3)
        self.floc_mix_power = pyunits.convert(self.floc_power_consumption * self.floc_mixers,
            to_units=pyunits.kilowatts)
        self.total_power = (self.rapid_mix_power + self.floc_mix_power) 
        return self.total_power / self.flow_in # kWh/m3

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2007
        tpec_tic = 'TPEC'
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.elect(),
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year, tpec_tic=tpec_tic)