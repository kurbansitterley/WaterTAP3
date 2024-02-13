from pyomo.environ import Var, Expression, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from watertap3.utils import cost_curves, financials
from watertap3.core.wt3_unit_pt import WT3UnitProcessPTData

## REFERENCE: ADD REFERENCE HERE

module_name = "basic_unit"


def cost_basic_unit(blk):
    basic_unit_name = blk.unit_model.config.unit_params["unit_process_name"]
    (
        blk.flow_basis,
        blk.cap_basis,
        blk.cap_exp,
        blk.energy_intensity,
        blk.basis_year,
        blk.base_kind,
    ) = cost_curves.basic_unit(basic_unit_name)
    blk.basis_currency = getattr(pyunits, f"MUSD_{blk.basis_year}")

    if blk.base_kind == "flow":
        blk.flow_scaling_base = Var(
            initialize=blk.flow_basis,
            bounds=(0, None),
            units=pyunits.m**3 / pyunits.hr,
            doc="Flow scaling base",
        )
        blk.flow_in = pyunits.convert(
            blk.unit_model.properties_in.flow_vol, to_units=pyunits.m**3 / pyunits.hr
        )

    elif blk.base_kind == "mass":
        blk.flow_scaling_base = Var(
            initialize=blk.flow_basis,
            bounds=(0, None),
            units=pyunits.kg / pyunits.hr,
            doc="Flow scaling base",
        )
        blk.flow_in = pyunits.convert(
            sum(
                blk.unit_model.properties_in.flow_mass_comp[j]
                for j in blk.unit_model.config.property_package.component_list
            ),
            to_units=pyunits.kg / pyunits.hr,
        )

    blk.capital_cost_base = Var(
        initialize=blk.cap_basis,
        bounds=(0, None),
        units=blk.basis_currency,
        doc=f"{basic_unit_name} (basic unit) capital cost base",
    )
    blk.capital_cost_exp = Var(
        initialize=blk.cap_exp,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc=f"{basic_unit_name} (basic unit) capital cost exponent",
    )
    blk.energy_intensity = Var(
        initialize=blk.energy_intensity,
        bounds=(0, None),
        units=pyunits.kWh / pyunits.m**3,
        doc=f"{basic_unit_name} (basic unit) energy intensity",
    )

    blk.fix_all_vars()

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)

    @blk.Constraint(doc=f"{basic_unit_name} (basic unit) capital cost equation")
    def capital_cost_constraint(b):
        flow_factor = b.flow_in / b.flow_scaling_base
        return b.capital_cost == pyunits.convert(
            b.capital_cost_base * flow_factor**b.capital_cost_exp,
            to_units=blk.costing_package.base_currency,
        )

    @blk.Expression(doc=f"{basic_unit_name} (basic unit) power required")
    def power_required(b):
        return pyunits.convert(
            b.energy_intensity * b.unit_model.properties_in.flow_vol,
            to_units=pyunits.kilowatt,
        )

    blk.costing_package.cost_flow(blk.power_required, "electricity")


@declare_process_block_class("UnitProcess")
class UnitProcessData(WT3UnitProcessPTData):
    def build(self):
        super().build()

    @property
    def default_costing_method(self):
        return cost_basic_unit

    def fixed_cap(self):
        """
        :param flow_in: Flow in to basic unit [m3/hr]
        :type flow_in: float
        """
        time = self.flowsheet().config.time.first()
        sys_cost_params = self.parent_block().costing_param
        flow_in_m3yr = pyunits.convert(
            self.flow_in, to_units=(pyunits.m**3 / pyunits.year)
        )

        if self.unit_process_name == "tramp_oil_tank":
            disposal_cost = 0.000114  # Kiran's disposal cost assumption $/m3
            self.costing.other_var_cost = (
                flow_in_m3yr
                * disposal_cost
                * sys_cost_params.plant_cap_utilization
                * 1e-6
            )

        if self.kind == "flow":
            flow_basis = self.basis * (pyunits.m**3 / pyunits.hour)
            flow_factor = self.flow_in / flow_basis
            basic_cap = self.cap_basis * flow_factor**self.cap_exp
            return basic_cap

        if self.kind == "mass":
            constituents = self.config.property_package.component_list
            mass_basis = self.basis * (pyunits.kg / pyunits.hour)
            mass_in = 0
            for constituent in constituents:
                mass_in += self.conc_mass_in[time, constituent]
            density = 0.6312 * mass_in + 997.86
            total_mass_in = density * self.flow_in
            mass_factor = total_mass_in / mass_basis
            basic_cap = self.cap_basis * mass_factor**self.cap_exp
            return basic_cap

    def elect(self):
        """
        Electricity intensity for basic units.

        :return: Electricity intensity [kWh/m3]
        """
        if (
            self.unit_process_name in ["mbr_nitrification", "mbr_denitrification"]
            and not self.case_specific
        ):
            # Electricity consumption for MBRs from:
            # "Assessing Location and Scale of Urban Nonpotable Water Reuse Systems for Life-Cycle Energy Consumption and Greenhouse Gas Emissions" Kavvada et al (2016)
            # Equation located in SI
            return 9.5 * self.flow_in**-0.3
        else:
            return self.elect_intensity

    def get_costing(self):
        """
        Initialize the unit in WaterTAP3.
        """
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(
            self.flow_vol_in[time], to_units=pyunits.m**3 / pyunits.hr
        )
        self.unit_process_name = self.unit_params["unit_process_name"]
        if "case_specific" in self.unit_params.keys():
            self.case_specific = self.unit_params["case_specific"]
            (
                self.basis,
                self.cap_basis,
                self.cap_exp,
                self.elect_intensity,
                self.basis_year,
                self.kind,
            ) = cost_curves.basic_unit(
                self.unit_process_name, case_specific=self.case_specific
            )
        else:
            self.case_specific = False
            (
                self.basis,
                self.cap_basis,
                self.cap_exp,
                self.elect_intensity,
                self.basis_year,
                self.kind,
            ) = cost_curves.basic_unit(self.unit_process_name)

        self.deltaP_outlet.unfix()
        self.deltaP_waste.unfix()
        self.pressure_out.fix(1)
        self.pressure_waste.fix(1)

        self.costing.fixed_cap_inv_unadjusted = Expression(
            expr=self.fixed_cap(), doc="Unadjusted fixed capital investment"
        )
        self.electricity = Expression(
            expr=self.elect(), doc="Electricity intensity [kWh/m3]"
        )
        financials.get_complete_costing(self.costing, basis_year=self.basis_year)
