import os
import pandas as pd
from pyomo.environ import Var, Constraint, Expression, Param, units as pyunits
from idaes.core import declare_process_block_class
from idaes.core.base.costing_base import (
    FlowsheetCostingBlockData,
    UnitModelCostingBlockData,
    register_idaes_currency_units,
)

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
case_study_TEA_basis_file = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/data/case_study_TEA_basis.csv"
)
industrial_electricity_costs_file = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/data/industrial_electricity_costs_2020.csv"
)
chem_costs_file = (
    os.path.abspath(os.path.join(__location__, os.pardir)) + "/data/chemical_costs.csv"
)
plant_costs_indices_file = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/data/plant_cost_indices.csv"
)


@declare_process_block_class("WT3Costing")
class WT3CostingData(FlowsheetCostingBlockData):
    def build_global_params(self):
        register_idaes_currency_units()
        case_study = self.parent_block().train["case_study"]
        self.base_currency = pyunits.USD_2021
        self.base_period = pyunits.year

        self.basis_data = pd.read_csv(case_study_TEA_basis_file, index_col="case_study").loc[case_study]
        self.basis_data.set_index("variable", inplace=True)
        self.basis_data = self.basis_data[["value"]]
        electricity_price_df = pd.read_csv(
            industrial_electricity_costs_file, index_col="location"
        )
        electricity_price_df.index = electricity_price_df.index.str.lower()
        # TODO: change these param names in csv
        basis_params = [
            "land_cost_percent",
            "working_capital_percent",
            "base_salary_per_fci",
            "maintenance_cost_percent",
            "laboratory_fees_percent",
            "insurance_and_taxes_percent",
            "employee_benefits_percent",
            "plant_life_yrs",
            "analysis_year",
            "debt_interest_rate",
            "plant_cap_utilization",
        ]

        self.location = (
            self.basis_data.loc["location_basis"].value
        )
        self.electricity_cost = Var(
            initialize=electricity_price_df.loc[self.location].iloc[0],
            units=self.base_currency / pyunits.kWh,
            doc=f"Electricity price for {self.location.title()}",
        )
        self.defined_flows["electricity"] = self.electricity_cost

        for bp in basis_params:
            p = float(
                self.basis_data.loc[bp].value
            )
            v = Var(
                initialize=p,
                doc=f"Cost basis parameter: {bp.replace('_', ' ').title()}",
            )
            self.add_component(bp, v)

        self.fix_all_vars()

@declare_process_block_class("WT3UnitCosting")
class WT3UnitCostingData(UnitModelCostingBlockData):
    pass