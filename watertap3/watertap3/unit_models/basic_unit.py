from pyomo.environ import Var, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from watertap3.utils import cost_curves
from watertap3.core.wt3_unit_pt import WT3UnitProcessPTData


module_name = "basic_unit"


def cost_basic_unit(blk):
    
    basic_unit_name = blk.unit_model.config.unit_params["unit_process_name"]
    (
        blk.flow_basis,
        blk.cap_basis,
        blk.cap_exp,
        blk._energy_intensity,
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
        initialize=blk._energy_intensity,
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


@declare_process_block_class("BasicUnit")
class UnitProcessData(WT3UnitProcessPTData):
    def build(self):
        super().build()

    @property
    def default_costing_method(self):
        return cost_basic_unit
