from pyomo.environ import Param, Var, Expression, Constraint, units as pyunits
from idaes.core.util.constants import Constants


def pumping_energy(blk):
    unit = blk.unit_model
    unit_params = unit.config.unit_params
    flow = unit.properties_in.flow_vol
    rho = unit.properties_in.dens_mass

    blk.pump_efficiency = Var(
        initialize=0.9,
        bounds=(0, 1),
        units=pyunits.dimensionless,
        doc="Pump efficiency",
    )
    blk.motor_efficiency = Var(
        initialize=0.9,
        bounds=(0, 1),
        units=pyunits.dimensionless,
        doc="Motor efficiency",
    )

    blk.lift_height = Var(
        initialize=100, bounds=(1, 1e5), units=pyunits.feet, doc="Pump lift height"
    )
    blk.pumping_power_required = Expression(
        expr=pyunits.convert(
            (flow * rho * Constants.acceleration_gravity * blk.lift_height)
            / (blk.pump_efficiency * blk.motor_efficiency),
            to_units=pyunits.kW,
        )
    )

    blk.fix_all_vars()

    blk.costing_package.cost_flow(blk.pumping_power_required, "electricity")
