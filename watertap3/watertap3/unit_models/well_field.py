from pyomo.environ import Var, Param, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from watertap3.core.wt3_unit_pt import WT3UnitProcessPTData

## REFERENCE: Derived from Voutchkov (2018) Table 4.6 and 4.7

module_name = "well_field"


def cost_well_field(blk):

    blk.basis_year = 2018
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")

    blk.capital_cost_base = Var(
        initialize=4731.6,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Well field capital cost basis",
    )
    blk.capital_cost_exp = Var(
        initialize=0.9196,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Well field capital cost exponent",
    )

    blk.piping_cost_basis = Var(
        initialize=35000,
        bounds=(0, None),
        units=blk.basis_currency / (pyunits.mile * pyunits.inch),
        doc="Piping cost basis",
    )

    if blk.unit_model.has_pump:
        blk.add_pumping_energy()
    blk.handle_costing_unit_params()

    blk.fix_all_vars()

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)
    flow_m3_hr = pyunits.convert(
        blk.unit_model.properties_in.flow_vol, to_units=pyunits.m**3 / pyunits.hr
    )

    blk.well_field_capital_cost = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Well field capital cost",
    )
    blk.piping_capital_cost = Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Piping capital cost",
    )

    @blk.Constraint(doc="Well field capital equation")
    def well_field_capital_cost_constraint(b):
        flow_dim = pyunits.convert(
            flow_m3_hr * (pyunits.hr / pyunits.m**3), to_units=pyunits.dimensionless
        )
        return (
            b.well_field_capital_cost
            == b.capital_cost_base * flow_dim**b.capital_cost_exp
        )

    @blk.Constraint(doc="Piping capital equation")
    def piping_capital_cost_constraint(b):
        return (
            b.piping_capital_cost
            == b.piping_cost_basis
            * b.unit_model.piping_distance
            * b.unit_model.piping_diameter
        )

    @blk.Constraint(doc="Unit total capital cost")
    def capital_cost_constraint(b):
        return b.capital_cost == pyunits.convert(
            b.well_field_capital_cost + b.piping_capital_cost,
            to_units=b.costing_package.base_currency,
        )


@declare_process_block_class("WellField")
class UnitProcessData(WT3UnitProcessPTData):
    def build(self):
        super().build()

        if (
            "pump" not in self.config.unit_params.keys()
            or not self.config.unit_params["pump"]
        ):
            self.has_pump = False
        else:
            self.has_pump = True

        self.piping_distance = Param(
            initialize=10, mutable=True, units=pyunits.mile, doc="Piping distance"
        )
        self.piping_diameter = Param(
            initialize=8, mutable=True, units=pyunits.inch, doc="Piping diameter"
        )
        self.handle_unit_params()

    @property
    def default_costing_method(self):
        return cost_well_field
