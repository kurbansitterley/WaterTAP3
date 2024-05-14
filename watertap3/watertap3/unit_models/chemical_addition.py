import os
import pandas as pd
from pyomo.environ import Var, Constraint, Param, Expression, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from watertap3.core.wt3_unit_pt import WT3UnitProcessPTData
from watertap3.core.util.pumping_energy import pumping_energy
from idaes.core.util.exceptions import ConfigurationError

## REFERENCE
## CAPITAL:
# Based on costs for SULFURIC ACID ADDITION 93% SOLUTION - FIGURE 5.5.11
# Cost Estimating Manual for Water Treatment Facilities (McGivney/Kawamura) (2008)
# DOI:10.1002/9780470260036
## ELECTRICITY:

module_name = "chemical_addition"
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
chemical_unit_cost_data_file = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/data/chemical_addition_cost_data.csv"
)


def cost_chemical_addition(blk):
    blk.basis_year = int(blk.unit_model.data.basis_year)
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")

    blk.capital_cost_base = Var(
        initialize=blk.unit_model.data.capital_cost_base,
        bounds=(0, None),
        units=blk.basis_currency,
        doc=f"{blk.unit_model.chemical.replace('_', ' ')} addition capital cost basis.",
    )
    blk.capital_cost_exp = Var(
        initialize=blk.unit_model.data.capital_cost_exp,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc=f"{blk.unit_model.chemical.replace('_', ' ')} addition capital cost exponent",
    )
    blk.add_pumping_energy(
        flow_vol=blk.unit_model.solution_volumetric_flow,
        rho=blk.unit_model.solution_density,
    )
    blk.handle_costing_unit_params()
    blk.fix_all_vars()

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TPEC")

    blk.cost_chemical_flow()

    @blk.Constraint(doc="Capital cost equation")
    def capital_cost_constraint(b):
        costing_flow_dim = pyunits.convert(
            b.unit_model.costing_flow * (b.unit_model.basis_units) ** -1,
            to_units=pyunits.dimensionless,
        )
        return b.capital_cost == pyunits.convert(
            (b.capital_cost_base * costing_flow_dim**b.capital_cost_exp)
            * b.unit_model.number_units
            * b.costing_package.TPEC,
            to_units=b.costing_package.base_currency,
        )


@declare_process_block_class("ChemicalAddition")
class UnitProcessData(WT3UnitProcessPTData):
    def build(self):
        super().build()

        if "chemical" not in self.config.unit_params.keys():
            raise ConfigurationError("Must provide a chemical name in unit_params.")
        else:
            chemical = self.config.unit_params["chemical"].lower()

        df = pd.read_csv(chemical_unit_cost_data_file, index_col="chemical_name")
        self.data = df.loc[chemical]
        self.costing_flow_basis = self.data.costing_flow_basis
        num, denom = [*self.data.basis_units.split("/")]
        self.basis_units = getattr(pyunits, num) / getattr(pyunits, denom)

        self.number_units = Param(
            initialize=2,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Number of units",
        )

        self.ratio_in_solution = Param(
            initialize=self.data.ratio_in_solution,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Ratio in solution",
        )

        self.solution_density = Param(
            initialize=self.data.solution_density,
            mutable=True,
            units=pyunits.kg / pyunits.m**3,
            doc="Chemical solution density",
        )

        self.dose = Param(
            initialize=0.01,
            mutable=True,
            units=pyunits.mg / pyunits.liter,
            doc="Chemical dose",
        )

        self.handle_unit_params()

        self.chemical = chemical

        if self.costing_flow_basis == "volume":
            self.get_vol_flow()
        if self.costing_flow_basis == "mass":
            self.get_mass_flow()

        # @self.Expression(doc=f"Mass flow of {self.chemical}")
        # def mass_flow_chemical(b):
        #     return pyunits.convert(b.properties_in.flow_vol, to_units=pyunits.kg/pyunits.s)

        @self.Expression(doc="Chemical volumetric flow used for pumping")
        def solution_volumetric_flow(b):
            return pyunits.convert(
                (b.properties_in.flow_vol * b.dose) / b.solution_density,
                to_units=pyunits.m**3 / pyunits.s,
            )

    def get_mass_flow(self):

        self.costing_flow = Var(
            initialize=10,
            bounds=(0, None),
            units=self.basis_units,
            doc="Mass flow rate used for capital cost.",
        )

        @self.Constraint(doc="Costing flow relationship")
        def costing_flow_constraint(b):
            cost_flow = b.properties_in.flow_vol * b.dose
            return b.costing_flow == pyunits.convert(
                cost_flow, to_units=self.basis_units
            )

    def get_vol_flow(self):

        self.costing_flow = Var(
            initialize=10,
            bounds=(0, None),
            units=self.basis_units,
            doc="Volmetric flow rate used for capital cost.",
        )

        @self.Constraint(doc="Costing flow relationship")
        def costing_flow_constraint(b):
            cost_flow = (
                b.properties_in.flow_vol * b.dose * b.ratio_in_solution
            ) / b.solution_density
            return b.costing_flow == pyunits.convert(
                cost_flow, to_units=self.basis_units
            )

    @property
    def default_costing_method(self):
        return cost_chemical_addition
