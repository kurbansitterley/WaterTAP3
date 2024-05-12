from pyomo.environ import Var, Param, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from watertap3.core.wt3_unit_pt import WT3UnitProcessPTData
from watertap3.core.wt3_unit_sido import WT3UnitProcessSIDOData

## REFERENCE
## CAPITAL:
# Kawamura and McGivney (2008)

module_name = "clarifier"


def cost_clarifier(blk):
    blk.basis_year = 2007
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")

    blk.capital_cost_base = Var(
        initialize=3470.6,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Clarifier capital cost base",
    )
    blk.capital_cost_exp = Var(
        initialize=0.6173,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Clarifier capital cost exponent",
    )

    blk.handle_costing_unit_params()
    blk.fix_all_vars()
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")

    @blk.Constraint(doc="Capital cost constraint")
    def capital_cost_constraint(b):
        floor_area_dim = pyunits.convert(
            b.unit_model.basin_floor_area * pyunits.ft**-2,
            to_units=pyunits.dimensionless,
        )
        return b.capital_cost == b.unit_model.number_clarifiers * pyunits.convert(
            b.capital_cost_base * floor_area_dim**b.capital_cost_exp,
            to_units=b.costing_package.base_currency,
        )


@declare_process_block_class("Clarifier")
class UnitProcessData(WT3UnitProcessSIDOData):
    def build(self):
        super().build()

        self.residence_time = Param(
            initialize=2, units=pyunits.hr, mutable=True, doc="Residence time"
        )

        self.basin_height = Param(
            initialize=10, units=pyunits.ft, mutable=True, doc="Basin height"
        )

        self.total_clarifier_volume = Var(
            initialize=1e3,
            units=pyunits.ft**3,
            bounds=(0, None),
            doc="Total clarifier volume",
        )

        self.number_clarifiers = Var(
            initialize=2,
            units=pyunits.dimensionless,
            bounds=(1, None),
            doc="Number of clarifiers",
        )

        self.basin_floor_area = Var(
            initialize=10000,
            units=pyunits.ft**2,
            bounds=(0, 35000),
            doc="Basin floor area",
        )

        self.handle_unit_params()

        @self.Expression(doc="Total clarifier floor area in ft2")
        def total_clarifier_floor_area(b):
            return pyunits.convert(
                b.total_clarifier_volume / b.basin_height, to_units=pyunits.ft**2
            )

        @self.Constraint(doc="Floor area for single clarifier")
        def basin_floor_area_constraint(b):
            return (
                b.basin_floor_area == b.total_clarifier_floor_area / b.number_clarifiers
            )

        @self.Constraint(doc="Total clarifier volume in ft3")
        def total_clarifier_volume_constraint(b):
            return b.total_clarifier_volume == pyunits.convert(
                b.properties_in.flow_vol * b.residence_time, to_units=pyunits.ft**3
            )

        @self.Constraint(doc="Number of clarifiers required")
        def number_clarifiers_constraint(b):
            return (
                b.number_clarifiers * b.basin_floor_area == b.total_clarifier_floor_area
            )

    @property
    def default_costing_method(self):
        return cost_clarifier
