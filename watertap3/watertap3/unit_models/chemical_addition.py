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
        flow_rate=blk.unit_model.solution_volumetric_flow,
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


# class UnitProcess(WT3UnitProcessPT):
#     def fixed_cap(self):
#         """
#         **"unit_params" are the unit parameters passed to the model from the input sheet as a Python dictionary.**

#         **EXAMPLE: {'dose': 10}**

#         :return: HCl addition fixed capital cost [$MM]
#         """
#         time = self.flowsheet().config.time.first()
#         self.flow_in = pyunits.convert(
#             self.flow_vol_in[time], to_units=pyunits.m**3 / pyunits.hr
#         )
#         self.number_of_units = 2
#         self.base_fixed_cap_cost = 900.97
#         self.cap_scaling_exp = 0.6179
#         chem_name = "Hydrochloric_Acid_(HCl)"
#         self.dose = Var(
#             initialize=1,
#             bounds=(0, None),
#             units=pyunits.kg / pyunits.m**3,
#             doc="Dose [kg/m3]",
#         )
#         self.dose.fix(0.010)
#         if "dose" in self.unit_params.keys():
#             self.dose.fix(self.unit_params["dose"] * 1e-3)
#         self.chem_dict = {chem_name: self.dose}
#         source_cost = (
#             self.base_fixed_cap_cost * self.solution_vol_flow() ** self.cap_scaling_exp
#         )
#         hcl_cap = (source_cost * self.tpec_tic * self.number_of_units) * 1e-6
#         return hcl_cap

#     def elect(self):
#         """
#         Electricity intensity.

#         :return: Electricity intensity [kWh/m3]
#         """
#         self.lift_height = 100 * pyunits.ft
#         self.pump_eff = 0.9 * pyunits.dimensionless
#         self.motor_eff = 0.9 * pyunits.dimensionless
#         soln_vol_flow = pyunits.convert(
#             self.solution_vol_flow(), to_units=(pyunits.gallon / pyunits.minute)
#         )
#         electricity = (
#             0.746
#             * soln_vol_flow
#             * self.lift_height
#             / (3960 * self.pump_eff * self.motor_eff)
#         ) / self.flow_in
#         return electricity

#     def solution_vol_flow(self):
#         """
#         Chemical solution flow in gal/day

#         :param solution_density: Solution density [kg/m3]
#         :type solution_density: float

#         :return: HCl solution flow [gal/day]
#         """
#         self.solution_density = 1490 * (pyunits.kg / pyunits.m**3)
#         chemical_rate = self.flow_in * self.dose
#         chemical_rate = pyunits.convert(
#             chemical_rate, to_units=(pyunits.kg / pyunits.day)
#         )
#         soln_vol_flow = chemical_rate / self.solution_density
#         soln_vol_flow = pyunits.convert(
#             soln_vol_flow, to_units=pyunits.gallon / pyunits.day
#         )
#         return soln_vol_flow

#     def get_costing(self):
#         """
#         Initialize the unit in WaterTAP3.
#         """
#         basis_year = 2007
#         tpec_tic = "TPEC"
#         self.costing.fixed_cap_inv_unadjusted = Expression(
#             expr=self.fixed_cap(), doc="Unadjusted fixed capital investment"
#         )
#         self.electricity = Expression(
#             expr=self.elect(), doc="Electricity intensity [kWh/m3]"
#         )
#         financials.get_complete_costing(
#             self.costing, basis_year=basis_year, tpec_tic=tpec_tic
#         )
