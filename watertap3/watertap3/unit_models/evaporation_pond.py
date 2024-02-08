from pyomo.environ import Expression, Constraint, \
    NonNegativeReals, Var, units as pyunits
from watertap3.utils import financials
from watertap3.core.wt3_unit_sido import WT3UnitProcess

## REFERENCE
## EVAPORATION RATE OF PURE WATER ESTIMATION:
# Turc (1961) PE in mm for day
# Jensen-Haise (1963) PE in mm per day
## CAPITAL:
## Costing for WT3 costing approach (the default approach) based on:
# Membrane Concentrate Disposal: Practices and Regulation (Second Edition) (2006) - Bureau of Reclamation
# Section 10 - Evaporation Pond Disposal
# usbr.gov/research/dwpr/reportpdfs/report123.pdf
##

module_name = 'evaporation_pond'

class UnitProcess(WT3UnitProcess):

    def fixed_cap(self):
        '''
        **"unit_params" are the unit parameters passed to the model from the input sheet as a Python dictionary.**
        Evaporation ponds can have many unit_params

        **EXAMPLE: {'approach': 'wt3', 'area': 3500, 'humidity': 0.75, 'wind_speed': 10}**

        :param unit_params: Input parameter dictionary from input sheet.
        :type unit_params: dict
        :param liner_thickness: Liner thickness [mil]
        :type liner_thickness: float
        :param land_cost: Cost of land for evaporation pond [$/acre]
        :type land_cost: float
        :param land_clearing_cost: Cost to clear land for evaporation pond [$/acre]
        :type land_clearing_cost: float
        :param dike_height: Height of dikes [ft]
        :type dike_height: float

        :return: Fixed capital cost for evaporation ponds [$MM]
        '''
        t = self.flowsheet().config.time.first()
        time = self.flowsheet().config.time
        
        self.tds_in = pyunits.convert(self.conc_mass_in[t, 'tds'],
            to_units=(pyunits.mg/pyunits.L))
        try:
            self.approach = self.unit_params['approach']
        except:
            self.approach = 'wt3'
        self.air_temp = Var(time,
            initialize=25,
            domain=NonNegativeReals,
            bounds=(0, 45),
            doc='Air Temp [C]')
        self.area = area = Var(time,
            initialize=1,
            domain=NonNegativeReals,
            units=pyunits.acres,
            doc='Evap. Pond Area [acres]')

        if 'area' in self.unit_params.keys():
            self.water_recovery.unfix()
            self.area.fix(self.unit_params['area'])

        try:
            self.liner_thickness = self.unit_params['liner_thickness']
            self.land_cost = self.unit_params['land_cost']
            self.land_clearing_cost = self.unit_params['land_clearing_cost']
            self.dike_height = self.unit_params['dike_height']
        except:
            self.liner_thickness = 50
            self.land_cost = 5000
            self.land_clearing_cost = 1000
            self.dike_height = 8

        self.evaporation_rate(t)
        self.evaporation_rate_regress(t)
        self.evap_rate = self.evap_rate_pure * 0.7
        self.flow_in = pyunits.convert(self.flow_vol_in[t],
            to_units=(pyunits.gallons/pyunits.minute))
        self.flow_waste = pyunits.convert(self.flow_vol_waste[t],
            to_units=(pyunits.gallons/pyunits.minute))
        self.flow_out = pyunits.convert(self.flow_vol_out[t],
            to_units=(pyunits.gallons/pyunits.minute))

        self.flow_constr = Constraint(expr=
            self.area[t] * self.evap_rate == self.flow_out)
        self.total_area = 1.2 * self.area[t] * (1 + 0.155 * self.dike_height / \
            (self.area[t] ** 0.5))
        self.cost_per_acre = 5406 + 465 * self.liner_thickness + 1.07 * self.land_cost + \
            0.931 * self.land_clearing_cost + 217.5 * self.dike_height
        if self.approach == 'zld':
            return 0.3 * area
        elif self.approach == 'lenntech':
            flow_in_m3_d = pyunits.convert(self.flow_in,
                to_units=(pyunits.m**3/pyunits.day))
            return 0.03099 * flow_in_m3_d ** 0.7613
        elif self.approach == 'wt3':
            return (self.cost_per_acre * self.total_area) * 1E-6


    def evaporation_rate(self, t):
        '''
        Calculation of evaporation rate [gpm/acre]

        :param unit_params: Input dictionary from input sheet.
        :type unit_params: dict
        :param t: Time indexing variable to use in Var()
        :type t: int
        :param evap_method: Evaporation rate method
        :type evap_method: str
        :param humidity: Humidity expressed as decimal for evaporation rate calculation
        :type humidity: float
        :param wind_speed: Wind speed for evaporation rate calculation [m/s]
        :type wind_speed: float
        :param air_temp: Air temperature for evaporation rate calculation [C]
        :type air_temp: float
        :param solar_rad: Incident solar radiation for evaporation rate calculation [mJ/m2]
        :type solar_rad: float


        :return:
        '''
        try:
            self.evap_method = self.unit_params['evap_method']
        except:
            self.evap_method = False
        try:
            self.humidity = self.unit_params['humidity']
            self.wind_speed = self.unit_params['wind_speed']
        except:
            self.humidity = 0.5
            self.wind_speed = 5
        if self.evap_method:

            try:
                self.air_temp.fix(self.unit_params['air_temp'])
                self.solar_rad = self.unit_params['solar_rad']
            except:
                self.air_temp.fix(20)
                self.solar_rad = 25
            if self.evap_method == 'turc':
                # Turc (1961) PE in mm for day
                self.evap_rate_pure = (0.313 * self.air_temp[t] * (self.solar_rad + 2.1) / \
                    (self.air_temp[t] + 15)) * (pyunits.millimeter / pyunits.day)
                self.evap_rate_pure = pyunits.convert(self.evap_rate_pure,
                to_units=(pyunits.gallons/pyunits.minute/pyunits.acre))
            if self.evap_method == 'jensen':
                # Jensen-Haise (1963) PE in mm per day
                self.evap_rate_pure = (0.41 * (0.025 * self.air_temp[t] + 0.078) * \
                    self.solar_rad) * (pyunits.millimeter/pyunits.day)
                self.evap_rate_pure = pyunits.convert(self.evap_rate_pure, 
                    to_units=(pyunits.gallons/pyunits.minute/pyunits.acre))
        else:
            # defaults to jensen
            self.air_temp.fix(25)
            self.solar_rad = 25  # average for 40deg latitude
            self.evap_rate_pure_mm_d = (0.41 * (0.025 * self.air_temp[t] + 0.078) * \
                self.solar_rad) * (pyunits.millimeter / pyunits.day)
            self.evap_rate_pure = pyunits.convert(self.evap_rate_pure_mm_d,
                to_units=(pyunits.gallons/pyunits.minute/pyunits.acre))
            self.evap_rate_m_yr = pyunits.convert(self.evap_rate_pure_mm_d,
                to_units=(pyunits.meter/pyunits.year))

        ## This costing model adapted from
        # Membrane Concentrate Disposal: Practices and Regulation (April 2006) - Bureau Land Management


    def evaporation_rate_regress(self, t):
        '''
        **NOTE: THIS FUNCTION IS NOT USED IN THE CURRENT RELEASE OF WaterTAP3**

        Calculates evaporation rate based on air temperature, TDS in, humidity, and wind speed.

        '''
        x0 = self.air_temp[t]
        x1 = self.tds_in
        x2 = self.humidity
        x3 = self.wind_speed
        ## CHANGED RATIO FUNCTION TO BE DEGREE=2 5/2/2021 -KAS
        ## THIS ISN'T USED CURRENTLY - NEEDS TO BE TESTED
        self.ratio = 0.0181657227 * (x0) + 4.38801e-05 * (x1) + 0.2504964875 * (x2) + 0.011328485 * (x3) - 0.0003463853 * (x0 ** 2) - 2.16888e-05 * (x0 * x1) - 0.0181098164 * (
                x0 * x2) + 0.0002098163 * (x0 * x3) + 8.654e-07 * (x1 ** 2) - 0.0004358946 * (x1 * x2) - 8.73918e-05 * (x1 * x3) - 0.0165224935 * (x2 ** 2) - 0.0174278724 * (
                             x2 * x3) + 0.0003850584 * (x3 ** 2) + 0.7943236298


    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        basis_year = 2007
        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=0,
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=basis_year)