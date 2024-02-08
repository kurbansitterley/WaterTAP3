# “Institute for the Design of Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE Framework) Copyright (c) 2019, by the software owners:
# The Regents of the University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University Research Corporation, et al.
# All rights reserved."
import idaes.logger as idaeslog
from idaes.core import declare_process_block_class
from pyomo.environ import Var, units as pyunits
from pyomo.network import Port
from .wt3_unit_base import WT3UnitProcessBaseData

__all__ = ["WT3UnitProcessPT"]


@declare_process_block_class("WT3UnitProcessPT")
class WT3UnitProcessPTData(WT3UnitProcessBaseData):
    """
    WaterTAP3 passthrough unit.
    Unit models using this class have no component removal and no waste stream.
    """
    def build(self):
        super().build()
        units_meta = self.config.property_package.get_metadata().get_derived_units

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True
        self.properties_in = prop_in = self.config.property_package.state_block_class(
            doc="Material properties of inlet stream", **tmp_dict
        )
        tmp_dict["defined_state"] = False
        self.properties_out = prop_out = self.config.property_package.state_block_class(
            doc="Material properties of outlet stream", **tmp_dict
        )

        self.deltaP = Var(
            initialize=0,
            units=units_meta("pressure"),
            doc="Pressure change between inlet and outlet",
        )

        self.deltaP.fix(0)

        @self.Constraint(doc="Pressure balance")
        def pressure_balance(b):
            return prop_in.pressure + b.deltaP == prop_out.pressure

        @self.Constraint(doc="Overall flow balance")
        def flow_balance(b):
            return prop_in.flow_vol == prop_out.flow_vol

        @self.Constraint(
            self.config.property_package.component_list, doc="Component mass balances"
        )
        def component_mass_balance(b, j):
            return prop_in.flow_mass_comp[j] == prop_out.flow_mass_comp[j]

        @self.Constraint(doc="Outlet temperature equation")
        def isothermal(b):
            return prop_in.temperature == prop_out.temperature
        
        self.inlet = Port(noruleinit=True, doc='Inlet Port')
        self.inlet.add(prop_in.flow_vol, 'flow_vol')
        self.inlet.add(prop_in.conc_mass_comp, 'conc_mass')
        self.inlet.add(prop_in.temperature, 'temperature')
        self.inlet.add(prop_in.pressure, 'pressure')

        self.outlet = Port(noruleinit=True, doc='Outlet Port')
        self.outlet.add(prop_out.flow_vol, 'flow_vol')
        self.outlet.add(prop_out.conc_mass_comp, 'conc_mass')
        self.outlet.add(prop_out.temperature, 'temperature')
        self.outlet.add(prop_out.pressure, 'pressure')