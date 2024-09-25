import math

from copy import deepcopy
from pyomo.environ import (
    Var,
    Param,
    exp,
    check_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

import idaes.logger as idaeslog
from idaes.core import declare_process_block_class
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.scaling import set_scaling_factor, get_scaling_factor

from watertap.costing.util import make_capital_cost_var, make_fixed_operating_cost_var
from watertap3.core.wt3_unit_pt import WT3UnitProcessPTData

## REFERENCE
## CAPITAL:
# Membrane Concentrate Disposal: Practices and Regulation, Chapter 10
# Mickley & Associates, 2006
# https://www.usbr.gov/research/dwpr/reportpdfs/report123.pdf
## ELECTRICITY:
# None

module_name = "evaporation_pond"


area_correction_factor_param_dict = {
    12: (2.6429, -0.202),
    8: (2.0512, -0.152),
    4: (1.5357, -0.092),
}

dike_cost_param_dict = {
    12: (29471, -0.317),
    8: (19000, -0.369),
    4: (84624, -0.431),
}

nominal_liner_cost_param_dict = {
    12: (1176, 0.2021),
    8: (15151, 0.1516),
    4: (20255, 0.0912),
}

fence_cost_param_dict = {
    12: (13505, -0.455),
    8: (11613, -0.423),
    4: (1024, -0.398),
}

road_cost_param_dict = {
    12: (2244, -0.454),
    8: (1943, -0.424),
    4: (1712.6, -0.399),
}


def cost_evaporation_pond(blk):
    blk.basis_year = 2001
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")
    evap_area_acre_dim = pyunits.convert(
        blk.unit_model.evaporative_area_acre * pyunits.acre**-1,
        to_units=pyunits.dimensionless,
    )

    total_area_required_acre = pyunits.convert(
        blk.unit_model.evaporation_pond_area * blk.unit_model.number_evaporation_ponds,
        to_units=pyunits.acre,
    )

    blk.liner_thickness_base = Param(
        initialize=60,
        units=pyunits.mil,
        doc="Basis for adjusting liner cost based on thickness",
    )

    blk.dike_capital_cost_base = Var(
        initialize=dike_cost_param_dict[blk.unit_model.dike_height][0],
        bounds=(0, None),
        units=blk.basis_currency / pyunits.acre,
        doc="Dike capital cost per acre base",
    )

    blk.dike_capital_cost_exp = Var(
        initialize=dike_cost_param_dict[blk.unit_model.dike_height][1],
        bounds=(None, 0),
        units=pyunits.dimensionless,
        doc="Dike capital cost per acre exponent",
    )

    blk.nominal_liner_capital_cost_base = Var(
        initialize=nominal_liner_cost_param_dict[blk.unit_model.dike_height][0],
        bounds=(0, None),
        units=blk.basis_currency / pyunits.acre,
        doc="Nominal liner capital cost per acre base",
    )

    blk.nominal_liner_capital_cost_exp = Var(
        initialize=nominal_liner_cost_param_dict[blk.unit_model.dike_height][1],
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Nominal liner capital cost per acre exponent",
    )

    blk.fence_capital_cost_base = Var(
        initialize=fence_cost_param_dict[blk.unit_model.dike_height][0],
        bounds=(0, None),
        units=blk.basis_currency / pyunits.acre,
        doc="Fence capital cost per acre base",
    )

    blk.fence_capital_cost_exp = Var(
        initialize=fence_cost_param_dict[blk.unit_model.dike_height][1],
        bounds=(None, 0),
        units=pyunits.dimensionless,
        doc="Fence capital cost per acre exponent",
    )

    blk.road_capital_cost_base = Var(
        initialize=road_cost_param_dict[blk.unit_model.dike_height][0],
        bounds=(0, None),
        units=blk.basis_currency / pyunits.acre,
        doc="Road capital cost per acre base",
    )

    blk.road_capital_cost_exp = Var(
        initialize=road_cost_param_dict[blk.unit_model.dike_height][1],
        bounds=(None, 0),
        units=pyunits.dimensionless,
        doc="Road capital cost per acre exponent",
    )

    blk.land_cost = Var(
        initialize=5000,
        bounds=(0, None),
        units=pyunits.USD_2001 / pyunits.acre,
        doc="Land cost per acre",
    )

    blk.land_clearing_cost = Var(
        initialize=5000,
        bounds=(0, None),
        units=pyunits.USD_2001 / pyunits.acre,
        doc="Land clearing cost per acre",
    )

    blk.liner_thickness = Var(
        initialize=60,
        bounds=(20, 120),
        units=pyunits.mil,
        doc="Liner thickness",
    )

    blk.liner_replacement_frequency = Var(
        initialize=20,
        bounds=(0, None),
        units=pyunits.year,
        doc="Liner replacement frequency",
    )

    blk.precipitate_handling_cost = Var(
        initialize=0,
        bounds=(0, None),
        units=blk.costing_package.base_currency / pyunits.kg,
        doc="Cost to handle precipitated solids",
    )

    blk.handle_costing_unit_params()
    blk.fix_all_vars()

    blk.dike_cost_per_acre = Var(
        initialize=5000,
        bounds=(0, None),
        units=blk.basis_currency / pyunits.acre,
        doc="Dike capital cost per acre",
    )

    blk.nominal_liner_cost_per_acre = Var(
        initialize=5000,
        bounds=(0, None),
        units=blk.basis_currency / pyunits.acre,
        doc="Nominal liner cost per acre",
    )

    blk.liner_cost_per_acre = Var(
        initialize=5000,
        bounds=(0, None),
        units=blk.basis_currency / pyunits.acre,
        doc="Liner cost per acre",
    )

    blk.fence_cost_per_acre = Var(
        initialize=5000,
        bounds=(0, None),
        units=blk.basis_currency / pyunits.acre,
        doc="Fence cost per acre",
    )

    blk.road_cost_per_acre = Var(
        initialize=5000,
        bounds=(0, None),
        units=pyunits.USD_2001 / pyunits.acre,
        doc="Road cost per acre",
    )

    ########

    blk.land_capital_cost = Var(
        initialize=500000,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Land capital cost",
    )

    blk.land_clearing_capital_cost = Var(
        initialize=500000,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Land clearing capital cost",
    )

    blk.dike_capital_cost = Var(
        initialize=500000,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Dike capital cost",
    )

    blk.liner_capital_cost = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Liner capital cost",
    )

    blk.fence_capital_cost = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Fence capital cost",
    )

    blk.road_capital_cost = Var(
        initialize=5000,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Road capital cost",
    )

    blk.precipitate_handling_operating_cost = Var(
        initialize=5000,
        bounds=(0, None),
        units=blk.costing_package.base_currency / pyunits.year,
        doc="Operating cost to handle precipitated solids",
    )

    blk.liner_replacement_operating_cost = Var(
        initialize=5000,
        bounds=(0, None),
        units=blk.basis_currency / pyunits.year,
        doc="Operating cost to replace liner",
    )

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)
    make_fixed_operating_cost_var(blk)

    capital_cost_expr = 0

    @blk.Constraint(doc="Capital cost for land")
    def land_capital_cost_constraint(b):
        return b.land_capital_cost == total_area_required_acre * b.land_cost

    capital_cost_expr += blk.land_capital_cost

    @blk.Constraint(doc="Capital cost for land clearing")
    def land_clearing_capital_cost_constraint(b):
        return (
            b.land_clearing_capital_cost
            == total_area_required_acre * b.land_clearing_cost
        )

    capital_cost_expr += blk.land_clearing_capital_cost

    @blk.Constraint(doc="Dike cost per acre")
    def dike_cost_per_acre_constraint(b):
        return (
            b.dike_cost_per_acre
            == b.dike_capital_cost_base * evap_area_acre_dim**b.dike_capital_cost_exp
        )

    @blk.Constraint(doc="Nominal liner cost per acre")
    def nominal_liner_cost_per_acre_constraint(b):
        return (
            b.nominal_liner_cost_per_acre
            == b.nominal_liner_capital_cost_base
            * evap_area_acre_dim**b.nominal_liner_capital_cost_exp
        )

    @blk.Constraint(doc="Liner cost per acre")
    def liner_cost_per_acre_constraint(b):
        return b.liner_cost_per_acre == b.nominal_liner_cost_per_acre * (
            b.liner_thickness / b.liner_thickness_base
        )

    @blk.Constraint(doc="Fence cost per acre")
    def fence_cost_per_acre_constraint(b):
        return (
            b.fence_cost_per_acre
            == b.fence_capital_cost_base
            * evap_area_acre_dim**b.fence_capital_cost_exp
        )

    @blk.Constraint(doc="Road cost per acre")
    def road_cost_per_acre_constraint(b):
        return (
            b.road_cost_per_acre
            == b.road_capital_cost_base * evap_area_acre_dim**b.road_capital_cost_exp
        )

    @blk.Constraint(doc="Dike capital cost")
    def dike_capital_cost_constraint(b):
        return b.dike_capital_cost == b.dike_cost_per_acre * total_area_required_acre

    capital_cost_expr += blk.dike_capital_cost

    @blk.Constraint(doc="Liner capital cost")
    def liner_capital_cost_constraint(b):
        return b.liner_capital_cost == b.liner_cost_per_acre * total_area_required_acre

    capital_cost_expr += blk.liner_capital_cost

    @blk.Constraint(doc="Fence capital cost")
    def fence_capital_cost_constraint(b):
        return b.fence_capital_cost == b.fence_cost_per_acre * total_area_required_acre

    capital_cost_expr += blk.fence_capital_cost

    @blk.Constraint(doc="Road capital cost")
    def road_capital_cost_constraint(b):
        return b.road_capital_cost == b.road_cost_per_acre * total_area_required_acre

    capital_cost_expr += blk.road_capital_cost

    @blk.Constraint(doc="Capital cost for evaporation pond")
    def capital_cost_constraint(b):
        return b.capital_cost == pyunits.convert(
            capital_cost_expr, to_units=blk.costing_package.base_currency
        )

    fixed_operating_cost_expr = 0

    @blk.Constraint(doc="Liner replacement cost")
    def liner_replacement_operating_cost_constraint(b):
        return b.liner_replacement_operating_cost == pyunits.convert(
            b.liner_capital_cost / b.liner_replacement_frequency,
            to_units=b.costing_package.base_currency / b.costing_package.base_period,
        )

    fixed_operating_cost_expr += blk.liner_replacement_operating_cost

    @blk.Constraint(doc="Solids handling operating cost")
    def precipitate_handling_operating_cost_constraint(b):
        return b.precipitate_handling_operating_cost == pyunits.convert(
            b.unit_model.mass_flow_precipitate * b.precipitate_handling_cost,
            to_units=b.costing_package.base_currency / b.costing_package.base_period,
        )

    fixed_operating_cost_expr += blk.precipitate_handling_operating_cost

    @blk.Constraint(doc="Fixed operating cost for evaporation pond")
    def fixed_operating_cost_constraint(b):
        return b.fixed_operating_cost == fixed_operating_cost_expr


@declare_process_block_class("EvaporationPond")
class UnitProcessData(WT3UnitProcessPTData):
    def build(self):
        super().build()

        if "dike_height" in self.config.unit_params.keys():
            self.dike_height = self.config.unit_params["dike_height"]  # ft
            if self.dike_height not in [4, 8, 12]:
                self.dike_height = 8
        else:
            self.dike_height = 8  # ft

        self.area_correction_factor_base = Param(
            initialize=area_correction_factor_param_dict[self.dike_height][0],
            mutable=True,
            units=pyunits.dimensionless,
            doc="Area correction factor base",
        )

        self.area_correction_factor_exp = Param(
            initialize=area_correction_factor_param_dict[self.dike_height][1],
            mutable=True,
            units=pyunits.dimensionless,
            doc="Area correction factor exponent",
        )

        self.daily_change_water_temperature = Param(
            initialize=1,
            mutable=True,
            units=pyunits.degK / pyunits.day,
            doc="Daily change in water temperature",
        )

        self.mw_ratio = Param(
            initialize=0.622,
            units=pyunits.dimensionless,
            doc="Ratio of molecular weight of water vapor to air",
        )

        self.arden_buck_coeff_a = Param(
            initialize=6.1121,
            units=pyunits.millibar,
            doc="Arden Buck equation for saturation vapor pressure, a coefficient",  # https://en.wikipedia.org/wiki/Arden_Buck_equation
        )

        self.arden_buck_coeff_b = Param(
            initialize=18.678,
            units=pyunits.dimensionless,
            doc="Arden Buck equation for saturation vapor pressure, b coefficient",  # https://en.wikipedia.org/wiki/Arden_Buck_equation
        )

        self.arden_buck_coeff_c = Param(
            initialize=257.14,
            units=pyunits.dimensionless,
            doc="Arden Buck equation for saturation vapor pressure, c coefficient",  # https://en.wikipedia.org/wiki/Arden_Buck_equation
        )

        self.arden_buck_coeff_d = Param(
            initialize=234.5,
            units=pyunits.dimensionless,
            doc="Arden Buck equation for saturation vapor pressure, d coefficient",  # https://en.wikipedia.org/wiki/Arden_Buck_equation
        )

        self.latent_heat_of_vaporization_intercept = Param(
            initialize=2500.78,
            units=pyunits.kilojoule / pyunits.kg,
            doc="Intercept of latent heat of vaporization equation",  # doi:10.1029/2008JD010174
        )

        self.latent_heat_of_vaporization_slope = Param(
            initialize=2.3601,
            units=pyunits.kilojoule / (pyunits.kg * pyunits.degK),
            doc="Slope of latent heat of vaporization equation",  # doi:10.1029/2008JD010174
        )

        self.heat_capacity_air = Param(
            initialize=1013,
            mutable=True,
            units=pyunits.joule / (pyunits.kg * pyunits.degK),
            doc="Specific heat capacity of dry air",  # doi:10.1029/2008JD010174
        )

        self.heat_capacity_water = Param(
            initialize=4186,
            mutable=True,
            units=pyunits.joule / (pyunits.kg * pyunits.degK),
            doc="Specific heat capacity of water",  # doi:10.1029/2008JD010174
        )

        ## solids precipitation rate is a function of TDS concentration [g / L]
        ## solids_precipitation_rate [ft / yr] = a1 * C**2 + a2 * C + intercept
        ## solids_precipitation_rate [ft / yr] = 4.12e-6 * C**2 + 1.92e-4 * C + 1.15e-3

        self.solids_precipitation_rate_a1 = Param(
            initialize=4.12e-6,
            mutable=True,
            units=pyunits.feet * pyunits.year**-1,
            doc="Solids precipitation rate a1 coefficient",
        )

        self.solids_precipitation_rate_a2 = Param(
            initialize=1.92e-4,
            mutable=True,
            units=pyunits.feet * pyunits.year**-1,
            doc="Solids precipitation rate a2 coefficient",
        )

        self.solids_precipitation_rate_intercept = Param(
            initialize=1.15e-3,
            mutable=True,
            units=pyunits.feet * pyunits.year**-1,
            doc="Solids precipitation rate intercept",
        )

        self.dens_solids = Param(
            initialize=2.16,
            mutable=True,
            units=pyunits.g / pyunits.cm**3,
            doc="Density of precipitated solids",
        )

        self.pressure_atm = Param(
            initialize=101325,
            mutable=True,
            units=pyunits.Pa,
            doc="Atmospheric pressure",
        )

        self.dens_vap = Param(
            initialize=1.293,
            mutable=True,
            units=pyunits.kg / pyunits.m**3,
            doc="Density of air",
        )

        self.evaporation_pond_depth = Param(
            initialize=18,
            mutable=True,
            units=pyunits.inches,
            doc="Depth of evaporation pond",
        )

        self.evaporation_rate_adjustment_factor = Param(
            initialize=0.7,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Factor to reduce evaporation rate as a results of higher salinity",
        )

        self.air_temperature = Param(
            initialize=300, mutable=True, units=pyunits.degK, doc="Air temperature"
        )

        self.air_temperature_C = self.air_temperature - 273.15 * pyunits.degK

        self.relative_humidity = Param(
            initialize=0.5,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Relative humidity",
        )

        self.net_solar_radiation = Param(
            initialize=150,
            mutable=True,
            units=pyunits.watt / pyunits.m**2,
            doc="Net incident solar radiation",  # net shortwave radiation - net longwave radiation
        )

        self.net_heat_flux = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.watt / pyunits.m**2,
            doc="Net heat flux out of system (water, soil)",
        )

        self.area_correction_factor = Var(
            initialize=1,
            bounds=(0.99, 3.2),
            units=pyunits.dimensionless,
            doc="Area correction factor",
        )

        self.saturation_vapor_pressure = Var(
            initialize=6,
            bounds=(0, None),
            units=pyunits.kPa,
            doc="Saturation vapor pressure at air temperature",
        )

        self.vapor_pressure = Var(
            initialize=4,
            bounds=(0, None),
            units=pyunits.kPa,
            doc="Vapor pressure at air temperature",
        )

        self.latent_heat_of_vaporization = Var(
            initialize=2500,
            bounds=(0, None),
            units=pyunits.kilojoule / pyunits.kg,
            doc="Latent heat of vaporization of water",
        )

        self.bowen_ratio = Var(
            initialize=0.3,
            bounds=(0, 2),
            units=pyunits.dimensionless,
            doc="Bowen ratio for BREB calculation of evaporation rate",
        )

        self.psychrometric_constant = Var(
            initialize=0.06,
            bounds=(0, None),
            units=pyunits.kPa * pyunits.degK**-1,
            doc="Psychrometric constant",
        )

        self.mass_flux_water_vapor = Var(
            initialize=1e-5,
            bounds=(0, 1e-3),
            units=pyunits.kg / (pyunits.m**2 * pyunits.s),
            doc="Mass flux of water vapor evaporated according using BREB method",
        )

        self.evaporation_rate = Var(
            initialize=0.03,
            bounds=(0, None),
            units=pyunits.m / pyunits.s,
            doc="Evaporation rate",
        )

        self.solids_precipitation_rate = Var(
            initialize=0.01,
            bounds=(0, None),
            units=pyunits.feet / pyunits.year,
            doc="Rate at which solids precipitate ",
        )

        self.total_evaporative_area_required = Var(
            initialize=100000,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Total evaporative area required",
        )

        self.evaporative_area_per_pond = Var(
            initialize=10000,
            bounds=(0, 405000),
            units=pyunits.m**2,
            doc="Evaporative area required per pond",
        )

        self.evaporation_pond_area = Var(
            initialize=10000,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Area of single evaporation pond",
        )

        self.number_evaporation_ponds = Var(
            initialize=1,
            bounds=(1, 25),
            units=pyunits.dimensionless,
            doc="Number of evaporation ponds",
        )

        self.handle_unit_params()

        @self.Expression()
        def evap_rate_mm_d(b):
            return pyunits.convert(
                b.evaporation_rate, to_units=pyunits.mm / pyunits.day
            )

        @self.Expression(doc="Total pond area in acres")
        def total_pond_area_acre(b):
            return pyunits.convert(
                b.evaporation_pond_area * b.number_evaporation_ponds,
                to_units=pyunits.acre,
            )

        @self.Expression(doc="Evaporative area per pond in acres")
        def evaporative_area_acre(b):
            return pyunits.convert(b.evaporative_area_per_pond, to_units=pyunits.acre)

        @self.Expression(doc="Individual pond area in acres")
        def pond_area_acre(b):
            return pyunits.convert(b.evaporation_pond_area, to_units=pyunits.acre)

        @self.Expression(doc="Net radiation for evaporation")
        def net_radiation(b):
            return b.net_solar_radiation - b.net_heat_flux

        @self.Expression(doc="Water temperature in deg C")
        def water_temperature_C(b):
            return b.properties_in.temperature - 273.15 * pyunits.degK

        @self.Expression(doc="Mass flow of precipitated solids")
        def mass_flow_precipitate(b):
            return pyunits.convert(
                b.number_evaporation_ponds
                * b.evaporation_pond_area
                * b.solids_precipitation_rate
                * b.dens_solids,
                to_units=pyunits.kg / pyunits.year,
            )

        @self.Expression(doc="Arden Buck fraction component")
        def arden_buck_exponential_term(b):
            temp_degC_dim = pyunits.convert(
                b.air_temperature_C * pyunits.degK**-1,
                to_units=pyunits.dimensionless,
            )  # dimensionless temperature in degC
            return (b.arden_buck_coeff_b - temp_degC_dim / b.arden_buck_coeff_d) * (
                temp_degC_dim / (b.arden_buck_coeff_c + temp_degC_dim)
            )

        @self.Constraint(doc="Solids precipitation rate")
        def solids_precipitation_rate_constraint(b):
            tds_in_dim = pyunits.convert(
                b.properties_in.conc_mass_comp["tds"] * pyunits.g**-1 * pyunits.L,
                to_units=pyunits.dimensionless,
            )
            return (
                b.solids_precipitation_rate
                == b.solids_precipitation_rate_a1 * tds_in_dim**2
                + b.solids_precipitation_rate_a2 * tds_in_dim
                + b.solids_precipitation_rate_intercept
            )

        @self.Constraint(doc="Area correction factor calculation")
        def area_correction_factor_constraint(b):
            evap_per_pond_acre_dim = pyunits.convert(
                b.evaporative_area_per_pond * pyunits.acre**-1,
                to_units=pyunits.dimensionless,
            )
            return (
                b.area_correction_factor
                == b.area_correction_factor_base
                * evap_per_pond_acre_dim**b.area_correction_factor_exp
            )

        @self.Constraint(doc="Water temperature")
        def water_temperature_constraint(b):
            return (
                b.properties_in.temperature
                == (1.04 * (b.air_temperature_C) + 0.22 * pyunits.degK)
                + 273.15 * pyunits.degK
            )

        @self.Constraint(doc="Net heat flux out of ecosystem")
        def net_flux_heat_constraint(b):
            return b.net_heat_flux == pyunits.convert(
                b.properties_in.dens_mass
                * b.heat_capacity_water
                * b.evaporation_pond_depth
                * b.daily_change_water_temperature,
                to_units=pyunits.watt / pyunits.m**2,
            )

        @self.Constraint(doc="Arden Buck equation for saturation vapor pressure")
        def arden_buck_constraint(b):
            return b.saturation_vapor_pressure == pyunits.convert(
                b.arden_buck_coeff_a * exp(b.arden_buck_exponential_term),
                to_units=pyunits.kPa,
            )

        @self.Constraint(doc="Latent heat of vaporization equation")
        def latent_heat_of_vaporization_constraint(b):
            return (
                b.latent_heat_of_vaporization
                == b.latent_heat_of_vaporization_intercept
                + b.latent_heat_of_vaporization_slope * b.air_temperature_C
            )

        @self.Constraint(doc="Psychrometric constant equation")
        def psychrometric_constant_constraint(b):
            return b.psychrometric_constant == pyunits.convert(
                (b.heat_capacity_air * b.pressure_atm)
                / (b.mw_ratio * b.latent_heat_of_vaporization),
                to_units=pyunits.kPa * pyunits.degK**-1,
            )

        @self.Constraint(doc="Vapor pressure equation")
        def vapor_pressure_constraint(b):
            return b.vapor_pressure == b.saturation_vapor_pressure * b.relative_humidity

        @self.Constraint(doc="Bowen ratio calculation")
        def bowen_ratio_constraint(b):
            return b.bowen_ratio == pyunits.convert(
                b.psychrometric_constant
                * (
                    (b.water_temperature_C - b.air_temperature_C)
                    / (b.saturation_vapor_pressure - b.vapor_pressure)
                ),
                to_units=pyunits.dimensionless,
            )

        @self.Constraint(doc="Mass flux water vapor")
        def mass_flux_water_vapor_constraint(b):
            return (
                b.mass_flux_water_vapor
                == b.evaporation_rate_adjustment_factor
                * pyunits.convert(
                    (b.net_radiation)
                    / (b.latent_heat_of_vaporization * (1 + b.bowen_ratio)),
                    to_units=pyunits.kg / (pyunits.m**2 * pyunits.s),
                )
            )

        @self.Constraint(doc="Evaporation rate")
        def evaporation_rate_constraint(b):
            return b.evaporation_rate == pyunits.convert(
                b.mass_flux_water_vapor / b.properties_in.dens_mass,
                to_units=pyunits.m / pyunits.s,
            )

        @self.Constraint(doc="Mass balance")
        def total_evaporative_area_required_constraint(b):
            return (
                b.mass_flux_water_vapor * b.total_evaporative_area_required
                == b.properties_in.flow_mass_comp["H2O"]
            )

        @self.Constraint(doc="Total evaporation pond area")
        def evaporative_area_per_pond_constraint(b):
            return (
                b.evaporative_area_per_pond * b.number_evaporation_ponds
                == b.total_evaporative_area_required
            )

        @self.Constraint(doc="Evaporation pond area")
        def evaporation_pond_area_constraint(b):
            return (
                b.evaporation_pond_area
                == b.evaporative_area_per_pond * b.area_correction_factor
            )

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for initialization routines

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        # return
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        if solver is None:
            opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        flags = self.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info("Initialization Step 1a Complete.")

        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # Set state_args from inlet state
        if state_args is None:
            self.state_args = state_args = {}
            state_dict = self.properties_in.define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        self.state_args_out = state_args_out = deepcopy(state_args)

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        cvc(self.properties_in.temperature, self.water_temperature_constraint)
        cvc(self.net_heat_flux, self.net_flux_heat_constraint)
        cvc(self.saturation_vapor_pressure, self.arden_buck_constraint)
        cvc(self.vapor_pressure, self.vapor_pressure_constraint)
        cvc(
            self.latent_heat_of_vaporization,
            self.latent_heat_of_vaporization_constraint,
        )
        cvc(self.psychrometric_constant, self.psychrometric_constant_constraint)
        cvc(self.bowen_ratio, self.bowen_ratio_constraint)
        cvc(self.mass_flux_water_vapor, self.mass_flux_water_vapor_constraint)

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)
                res = opt.solve(self, tee=slc.tee)
        n_ponds = math.floor(value(self.number_evaporation_ponds))
        self.number_evaporation_ponds.set_value(n_ponds)
        res = opt.solve(self, tee=slc.tee)

        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Release Inlet state
        self.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize.")

        self.initialized = True

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if get_scaling_factor(self.psychrometric_constant) is None:
            set_scaling_factor(self.psychrometric_constant, 1e2)

        if get_scaling_factor(self.bowen_ratio) is None:
            set_scaling_factor(self.bowen_ratio, 1e2)

        if get_scaling_factor(self.mass_flux_water_vapor) is None:
            set_scaling_factor(self.mass_flux_water_vapor, 1e4)

    @property
    def default_costing_method(self):
        return cost_evaporation_pond
