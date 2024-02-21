from pyomo.environ import Var, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from watertap3.core.wt3_unit_sido import WT3UnitProcessSIDOData

## REFERENCE
## CAPITAL:
## ELECTRICITY:
# Equation derived from:
# L.H. Shaffer and M.S. Mintz, Electrodialysis, in Principles of Desalination,
# K.S. Spiegler (ed.), Academic Press, New York, pp. 200-289 (1966).

module_name = "electrodialysis_reversal"


def cost_electrodialysis_reversal(blk):

    blk.basis_year = 2016
    blk.basis_currency = getattr(pyunits, f"MUSD_{blk.basis_year}")

    blk.flow_scaling_base = Var(
        initialize=946,
        bounds=(0, None),
        units=pyunits.m**3 / pyunits.hr,
        doc="Flow scaling base",
    )

    blk.capital_cost_slope = Var(
        initialize=31,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Capital cost slope",
    )

    blk.energy_intensity_intercept = Var(
        initialize=0.2534,
        bounds=(0, None),
        units=pyunits.kWh / pyunits.m**3,
        doc="Energy intensity equation intercept",
    )

    blk.energy_intensity_slope = Var(
        initialize=5.149e-4,
        bounds=(0, None),
        units=(pyunits.kWh * pyunits.liter) / (pyunits.m**3 * pyunits.mg),
        doc="Energy intensity equation slope",
    )

    blk.handle_costing_unit_params()
    blk.fix_all_vars()

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)

    blk.energy_intensity = Var(
        initialize=1,
        bounds=(0, None),
        units=pyunits.kWh / pyunits.m**3,
        doc="Energy intensity",
    )

    @blk.Expression(doc="Power required")
    def power_required(b):
        return pyunits.convert(
            b.energy_intensity * b.unit_model.properties_in.flow_vol,
            to_units=pyunits.kW,
        )

    @blk.Constraint(doc="Capital cost equation")
    def capital_cost_constraint(b):
        flow_in = pyunits.convert(
            b.unit_model.properties_in.flow_vol, to_units=pyunits.m**3 / pyunits.hr
        )
        return b.capital_cost == pyunits.convert(
            b.capital_cost_slope * flow_in / b.flow_scaling_base,
            to_units=b.costing_package.base_currency,
        )

    @blk.Constraint(doc="Energy intensiy equation")
    def energy_intensity_constraint(b):
        tds_in = pyunits.convert(
            b.unit_model.properties_in.conc_mass_comp["tds"],
            to_units=pyunits.mg / pyunits.liter,
        )
        return (
            b.energy_intensity
            == b.energy_intensity_slope * tds_in + b.energy_intensity_intercept
        )


@declare_process_block_class("ElectrodialysisReversal")
class UnitProcessData(WT3UnitProcessSIDOData):
    def build(self):
        super().build()

    @property
    def default_costing_method(self):
        return cost_electrodialysis_reversal


# class UnitProcess(WT3UnitProcess):
#     def fixed_cap(self):
#         time = self.flowsheet().config.time.first()
#         self.flow_in = pyunits.convert(
#             self.flow_vol_in[time], to_units=pyunits.m**3 / pyunits.hr
#         )

#         ed_cap = 31 * self.flow_in / 946  # $MM
#         return ed_cap

#     def elect(self):

#         # EDR electricity intensity updated June 23, 2021
#         # Correspondence with Cameron McKay on June 22, 2021
#         # This regression was done in excel EDR_electricity.xlsx
#         # Figure is from Membrane Technology and Applications, by Richard W. Baker and is citing the following source for this figure:
#         # L.H. Shaffer and M.S. Mintz, Electrodialysis, in Principles of Desalination, K.S. Spiegler (ed.), Academic Press, New York, pp. 200-289 (1966).
#         time = self.flowsheet().config.time.first()
#         self.tds_in = pyunits.convert(
#             self.conc_mass_in[time, "tds"], to_units=(pyunits.mg / pyunits.L)
#         )
#         ed_elect = 5.149e-4 * self.tds_in + 0.2534
#         return ed_elect

#     def get_costing(self):
#         """
#         Initialize the unit in WaterTAP3.
#         """
#         basis_year = 2016
#         self.costing.fixed_cap_inv_unadjusted = Expression(
#             expr=self.fixed_cap(), doc="Unadjusted fixed capital investment"
#         )
#         self.electricity = Expression(
#             expr=self.elect(), doc="Electricity intensity [kWh/m3]"
#         )
#         financials.get_complete_costing(self.costing, basis_year=basis_year)
