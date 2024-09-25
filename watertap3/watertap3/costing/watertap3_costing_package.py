import os
import pandas as pd
from pyomo.environ import Block, Var, Constraint, Expression, Param, units as pyunits
from idaes.core import declare_process_block_class
from idaes.core.util.constants import Constants
from idaes.core.base.costing_base import (
    UnitModelCostingBlockData,
    register_idaes_currency_units,
)
from idaes.core.util.scaling import set_scaling_factor, get_scaling_factor
from watertap.costing.costing_base import WaterTAPCostingBlockData
from watertap3.utils.financials import get_ind_table

__author__ = "Kurban Sitterley"
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
case_study_TEA_basis_file = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/data/case_study_TEA_basis.csv"
)
industrial_electricity_costs_file = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/data/industrial_electricity_costs_2020.csv"
)
material_costs_file = (
    os.path.abspath(os.path.join(__location__, os.pardir)) + "/data/chemical_costs.csv"
)
plant_costs_indices_file = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/data/plant_cost_indices.csv"
)
last_year_for_cost_indicies = 2050


@declare_process_block_class("WT3Costing")
class WT3CostingData(WaterTAPCostingBlockData):
    def build_global_params(self):
        register_idaes_currency_units()
        pyunits.load_definitions_from_strings(
            ["USD_2022 = 500/816.0 * USD_CE500", "USD_2023 = 500/797.9 * USD_CE500"]
        )
        case_study = self.parent_block().train["case_study"]
        self.base_period = pyunits.year

        self.basis_data = pd.read_csv(
            case_study_TEA_basis_file, index_col="case_study"
        ).loc[case_study]
        self.material_costs_data = pd.read_csv(
            material_costs_file, index_col="material"
        )
        self.basis_data.set_index("variable", inplace=True)
        self.basis_data = self.basis_data[["value"]]
        electricity_price_df = pd.read_csv(
            industrial_electricity_costs_file, index_col="location"
        )
        electricity_price_df.index = electricity_price_df.index.str.lower()
        self.basis_year = int(self.basis_data.loc["basis_year"].value)
        self.base_currency = getattr(pyunits, f"USD_{self.basis_year}")
        basis_params = [
            "land_cost_factor",  # land_cost_percent, land_cost_percent_FCI
            "working_capital_factor",  # working_capital_percent, working_cap_percent_FCI
            "salary_cost_factor",  # base_salary_per_fci, salaries_percent_FCI
            "maintenance_cost_factor",  # maintenance_cost_percent, maintenance_costs_percent_FCI
            "laboratory_cost_factor",  # laboratory_fees_percent, lab_fees_percent_FCI
            "insurance_taxes_cost_factor",  # insurance_and_taxes_percent, insurance_taxes_percent_FCI
            "benefits_cost_factor",  # employee_benefits_percent, benefit_percent_of_salary
            "plant_lifetime",
            "debt_interest_rate",
            "utilization_factor",
        ]

        self.location = self.basis_data.loc["location_basis"].value
        self.electricity_cost = Var(
            initialize=electricity_price_df.loc[self.location].iloc[0],
            units=self.base_currency / pyunits.kWh,
            doc=f"Electricity price for {self.location.title()}",
        )
        self.defined_flows["electricity"] = self.electricity_cost

        for bp in basis_params:
            p = float(self.basis_data.loc[bp].value)
            v = Var(
                initialize=p,
                doc=f"Cost basis parameter: {bp.replace('_', ' ').title()}",
            )
            self.add_component(bp, v)

        self.factor_total_investment = Var(
            initialize=1.0,
            doc="Total investment factor [investment cost/equipment cost]",
            units=pyunits.dimensionless,
        )

        self.capital_recovery_factor = Expression(
            expr=(
                self.debt_interest_rate
                * (1 + self.debt_interest_rate) ** self.plant_lifetime
            )
            / (((1 + self.debt_interest_rate) ** self.plant_lifetime) - 1)
        )

        self.TPEC = Var(
            initialize=3.4,
            units=pyunits.dimensionless,
            doc=f"Total purchased equipment cost factor (TPEC)",
        )
        self.TIC = Var(
            initialize=1.65,
            units=pyunits.dimensionless,
            doc=f"Total installed cost factor (TIC)",
        )
        self.fix_all_vars()

    def build_process_costs(self):
        # add total_captial_cost and total_operating_cost
        self._build_common_process_costs()
        # self.material_costs_data = pd.read_csv(material_costs_file, index_col="material")

        self.land_cost = Expression(
            expr=self.aggregate_capital_cost * self.land_cost_factor
        )

        self.working_capital_cost = Expression(
            expr=self.aggregate_capital_cost * self.working_capital_factor
        )

        self.salary_cost = Expression(
            expr=self.aggregate_capital_cost * self.salary_cost_factor
        )

        self.benefit_cost = Expression(
            expr=self.salary_cost * self.benefits_cost_factor
        )

        self.maintenance_cost = Expression(
            expr=self.aggregate_capital_cost * self.maintenance_cost_factor
        )

        self.lab_cost = Expression(
            expr=self.aggregate_capital_cost * self.laboratory_cost_factor
        )

        self.insurance_taxes_cost = Expression(
            expr=self.aggregate_capital_cost * self.insurance_taxes_cost_factor
        )

        self.total_fixed_operating_cost = Expression(
            expr=self.aggregate_fixed_operating_cost
            + self.salary_cost
            + self.benefit_cost
            + self.maintenance_cost
            + self.lab_cost
            + self.insurance_taxes_cost,
            doc="Total fixed operating costs",
        )

        self.total_variable_operating_cost = Expression(
            expr=(
                self.aggregate_variable_operating_cost
                + sum(self.aggregate_flow_costs[f] for f in self.used_flows)
                * self.utilization_factor
            )
            if self.used_flows
            else self.aggregate_variable_operating_cost,
            doc="Total variable operating cost of process per operating period",
        )

        self.total_capital_cost_constraint = Constraint(
            expr=self.total_capital_cost
            == self.factor_total_investment * self.aggregate_capital_cost
            + self.land_cost
            + self.working_capital_cost
        )

        self.total_operating_cost_constraint = Constraint(
            expr=self.total_operating_cost
            == (self.total_fixed_operating_cost + self.total_variable_operating_cost),
            doc="Total operating cost of process per operating period",
        )

        self.total_annualized_cost = Expression(
            expr=(
                self.total_capital_cost * self.capital_recovery_factor
                + self.total_operating_cost
            ),
            doc="Total annualized cost of operation",
        )

    def add_material_cost_param_block(self, material=None):
        if material is None:
            raise ValueError("Must provide chemical name to be costed.")
        mat_ser = self.material_costs_data.loc[material]
        n = getattr(pyunits, f"USD_{mat_ser.price_year}")
        d = getattr(pyunits, mat_ser.price_units.split("/")[-1])

        def build_rule(blk):
            blk.cost = Param(
                initialize=mat_ser.price,
                mutable=True,
                units=n / d,
                doc=f"{material.title().replace('_', ' ')} cost",
            )
            blk.purity = Param(
                initialize=mat_ser.purity,
                mutable=True,
                units=pyunits.dimensionless,
                doc=f"{material.title().replace('_', ' ')} purity",
            )
            blk.parent_block().defined_flows[material] = blk.cost
            blk.parent_block().register_flow_type(material, blk.cost / blk.purity)

        mat_block = Block(rule=build_rule)
        self.add_component(material, mat_block)


@declare_process_block_class("WT3UnitCosting")
class WT3UnitCostingData(UnitModelCostingBlockData):
    def build(self):
        super().build()

    def cost_chemical_flow(self, **kwargs):
        if "chemical" in kwargs.keys():
            kwargs["material"] = kwargs["chemical"]
        self.cost_material_flow(**kwargs)

    def cost_material_flow(
        self, material=None, flow_rate=None, name=None, dose=None, flow_expr=None
    ):
        if material is None:
            if hasattr(self.unit_model, "chemical"):
                material = self.unit_model.chemical
            elif hasattr(self.unit_model, "material"):
                material = self.unit_model.material
            else:
                raise ValueError("Must include a material name.")
        if material not in self.costing_package.flow_types:
            self.costing_package.add_material_cost_param_block(material=material)

        if flow_expr is None:
            if dose is None:
                dose = self.unit_model.dose
            if flow_rate is None:
                flow_rate = self.unit_model.properties_in.flow_vol
            flow_expr = Expression(
                expr=pyunits.convert(
                    flow_rate * dose,
                    to_units=pyunits.kg / self.costing_package.base_period,
                )
            )
        else:
            flow_expr = Expression(
                expr=pyunits.convert(
                    flow_expr, to_units=pyunits.kg / self.costing_package.base_period
                )
            )

        if name is None:
            name = f"{material}_flow"

        self.add_component(name, flow_expr)

        self.costing_package.cost_flow(flow_expr, material)

    def handle_costing_unit_params(self):
        for k, v in self.unit_model.config.unit_params.items():
            if hasattr(self, k):
                p = getattr(self, k)
                if isinstance(p, Var):
                    p.fix(v)
                elif isinstance(p, Param):
                    p.set_value(v)
                else:
                    raise ValueError(
                        f"{k} in unit_params for {self.name} is for a {type(p)}"
                        "but must be for a Var or Param."
                        f"Remove {k} from unit_params on the input sheet."
                    )

    def add_pumping_energy(
        self, flow_vol=None, flow_mass=None, rho=None, lift_height=100
    ):

        unit = self.unit_model
        if flow_vol is None:
            flow_vol = unit.properties_in.flow_vol
        if rho is None:
            rho = unit.properties_in.dens_mass
        if flow_mass is None:
            flow_mass = flow_vol * rho

        self.pump_efficiency = Var(
            initialize=0.9,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Pump efficiency",
        )
        self.motor_efficiency = Var(
            initialize=0.9,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Motor efficiency",
        )

        self.lift_height = Var(
            initialize=lift_height,
            bounds=(1, 1e5),
            units=pyunits.feet,
            doc="Pump lift height",
        )
        self.power_required = Expression(
            expr=pyunits.convert(
                (flow_mass * Constants.acceleration_gravity * self.lift_height)
                / (self.pump_efficiency * self.motor_efficiency),
                to_units=pyunits.kW,
            )
        )

        self.costing_package.cost_flow(self.power_required, "electricity")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # for v in self.component_objects(Var, descend_into=False):
        #     set_scaling_factor(v, 1e-5)

    def get_cost_indices(self):
        # if self.basis_year is None:
        #     return

        df = get_ind_table(self.costing_package.basis_year)
        self.capital_factor = Param(
            initialize=df.loc[self.basis_year].capital_factor,
            mutable=True,
            doc="Capital factor",
        )
        self.chemical_factor = Param(
            initialize=df.loc[self.basis_year].chemical_factor,
            mutable=True,
            doc="Chemical factor",
        )
        self.labor_factor = Param(
            initialize=df.loc[self.basis_year].labor_factor,
            mutable=True,
            doc="Labor factor",
        )
