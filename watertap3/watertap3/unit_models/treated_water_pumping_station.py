from pyomo.environ import Var, Constraint, Param, Expression, units as pyunits
from idaes.core import declare_process_block_class
from watertap.costing.util import make_capital_cost_var
from watertap3.core.wt3_unit_pt import WT3UnitProcessPTData
from watertap3.core.util.pumping_energy import pumping_energy

## REFERENCE: McGivney & Kawamura (2008) - Figure 5.5.35b

module_name = "municipal_drinking"


def cost_municipal_drinking(blk):
    blk.basis_year = 2007
    blk.basis_currency = getattr(pyunits, f"USD_{blk.basis_year}")

    blk.capital_cost_base = Var(
        initialize=45314,
        bounds=(0, None),
        units=blk.basis_currency,
        doc="Municipal drinking capital cost basis",
    )

    blk.capital_cost_exp = Var(
        initialize=0.84,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Municipal drinking capital cost exponent",
    )

    blk.add_pumping_energy()
    blk.lift_height.fix(300)
    blk.handle_costing_unit_params()

    blk.fix_all_vars()

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TPEC")
    flow_in_mgd = pyunits.convert(
        blk.unit_model.properties_in.flow_vol, to_units=pyunits.Mgallons / pyunits.day
    )

    @blk.Constraint(doc="Capital cost equation")
    def capital_cost_constraint(b):
        flow_in_dim = pyunits.convert(
            flow_in_mgd * pyunits.day * pyunits.Mgallons**-1,
            to_units=pyunits.dimensionless,
        )
        return b.capital_cost == pyunits.convert(
            b.capital_cost_base * flow_in_dim**b.capital_cost_exp,
            to_units=blk.costing_package.base_currency,
        )


@declare_process_block_class("TreatedWaterPumpingStation")
class UnitProcessData(WT3UnitProcessPTData):
    def build(self):
        super().build()

    @property
    def default_costing_method(self):
        return cost_municipal_drinking

    # class UnitProcess(WT3UnitProcessPT):

    def fixed_cap(self):
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(
            self.flow_vol_in[time], to_units=pyunits.Mgallons / pyunits.day
        )
        self.base_fixed_cap_cost = 0.0403
        self.cap_scaling_exp = 0.8657
        muni_cap = (
            self.base_fixed_cap_cost * self.flow_in**self.cap_scaling_exp
        ) * self.tpec_tic
        return muni_cap

    def elect(self):
        self.lift_height = 300 * pyunits.ft
        self.pump_eff = 0.9 * pyunits.dimensionless
        self.motor_eff = 0.9 * pyunits.dimensionless
        flow_in_gpm = pyunits.convert(
            self.flow_in, to_units=pyunits.gallons / pyunits.minute
        )
        flow_in_m3hr = pyunits.convert(
            self.flow_in, to_units=pyunits.m**3 / pyunits.hour
        )
        electricity = (
            0.746
            * flow_in_gpm
            * self.lift_height
            / (3960 * self.pump_eff * self.motor_eff)
        ) / flow_in_m3hr
        return electricity

    def get_costing(self):
        """
        Initialize the unit in WaterTAP3.
        """
        basis_year = 2020
        tpec_tic = "TPEC"
        self.costing.fixed_cap_inv_unadjusted = Expression(
            expr=self.fixed_cap(), doc="Unadjusted fixed capital investment"
        )
        self.electricity = Expression(
            expr=self.elect(), doc="Electricity intensity [kWh/m3]"
        )
        financials.get_complete_costing(
            self.costing, basis_year=basis_year, tpec_tic=tpec_tic
        )
