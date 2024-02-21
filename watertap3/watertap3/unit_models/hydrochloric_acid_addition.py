from pyomo.environ import Var, Constraint, Param, Expression, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from watertap3.core.wt3_unit_pt import WT3UnitProcessPTData
from watertap3.core.util.pumping_energy import pumping_energy

## REFERENCE
## CAPITAL:
# Based on costs for SULFURIC ACID ADDITION 93% SOLUTION - FIGURE 5.5.11
# Cost Estimating Manual for Water Treatment Facilities (McGivney/Kawamura) (2008)
# DOI:10.1002/9780470260036
## ELECTRICITY:

module_name = "hydrochloric_acid_addition"


def cost_hcl_addition(blk):
    blk.basis_year = 2007
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")

    blk.capital_cost_base = Var(
        initialize=900.97,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="HCl addition capital cost basis",
    )
    blk.capital_cost_exp = Var(
        initialize=0.6179,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="HCl addition capital cost exponent",
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
        sol_flow_dim = pyunits.convert(
            b.unit_model.solution_volumetric_flow * pyunits.day * pyunits.gallons**-1,
            to_units=pyunits.dimensionless,
        )
        return b.capital_cost == pyunits.convert(
            (b.capital_cost_base * sol_flow_dim**b.capital_cost_exp)
            * b.unit_model.number_units
            * b.costing_package.TPEC,
            to_units=b.costing_package.base_currency,
        )


@declare_process_block_class("HClAddition")
class UnitProcessData(WT3UnitProcessPTData):
    def build(self):
        super().build()

        if "chemical" in self.config.unit_params.keys():
            self.chemical = self.config.unit_params["chemical"]
        else:
            self.chemical = "hydrochloric_acid"

        self.number_units = Param(
            initialize=2,
            units=pyunits.dimensionless,
            doc="Number of units",
        )

        self.solution_density = Param(
            initialize=1490,
            units=pyunits.kg / pyunits.m**3,
            doc="Solution density",
        )

        self.dose = Param(
            initialize=0.01,
            units=pyunits.mg / pyunits.liter,
            doc="Chemical dose",
        )

        self.handle_unit_params()

        @self.Expression(doc="HCl volumetric flow")
        def solution_volumetric_flow(b):
            return pyunits.convert(
                (b.properties_in.flow_vol * b.dose) / b.solution_density,
                to_units=pyunits.gallon / pyunits.day,
            )

    @property
    def default_costing_method(self):
        return cost_hcl_addition
    





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
