# “Institute for the Design of Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE Framework) Copyright (c) 2019, by the software owners:
# The Regents of the University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University Research Corporation, et al.
# All rights reserved."
import idaes.logger as idaeslog
from idaes.core import declare_process_block_class
from pyomo.environ import Var, units as pyunits
from pyomo.network import Port
from .wt3_unit_base import WT3UnitProcessBaseData

__all__ = ["WT3UnitProcessSIDO"]


@declare_process_block_class("WT3UnitProcessSIDO")
class WT3UnitProcessSIDOData(WT3UnitProcessBaseData):
    """
    WaterTAP3 unit model with a single inlet that is split between and outlet and waste stream.
    Unit models using this class can have component removal and water recovery.
    Any component removal goes to the waste stream.
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
        self.properties_waste = (
            prop_waste
        ) = self.config.property_package.state_block_class(
            doc="Material properties of waste stream", **tmp_dict
        )

        self.deltaP_outlet = Var(
            initialize=1e-6,
            units=units_meta("pressure"),
            doc="Pressure change between inlet and outlet",
        )

        self.deltaP_outlet.fix(0)

        self.deltaP_waste = Var(
            initialize=1e-6,
            units=units_meta("pressure"),
            doc="Pressure change between inlet and waste",
        )
        self.deltaP_waste.fix(0)

        ## WATER RECOVERY & REMOVAL FRACTION
        self.water_recovery = Var(
            initialize=0.8,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Water recovery fraction",
        )
        self.removal_fraction = Var(
            self.config.property_package.component_list,
            initialize=0.01,
            units=pyunits.dimensionless,
            doc="Component removal fraction",
        )

        @self.Constraint(doc="Outlet pressure equation")
        def pressure_balance_outlet(b):
            return prop_in.pressure + b.deltaP_outlet == prop_out.pressure

        @self.Constraint(doc="Waste pressure equation")
        def pressure_balance_waste(b):
            return prop_in.pressure + b.deltaP_waste == prop_waste.pressure

        @self.Constraint(doc="Water recovery equation")
        def recovery_equation(b):
            return b.water_recovery * prop_in.flow_vol == prop_out.flow_vol

        @self.Constraint(doc="Overall flow balance")
        def flow_balance(b):
            return prop_in.flow_vol == prop_out.flow_vol + prop_waste.flow_vol

        @self.Constraint(
            self.config.property_package.component_list,
            doc="Component removal equation",
        )
        def component_removal_equation(b, j):
            return (
                b.removal_fraction[j] * prop_in.flow_mass_comp[j]
                == prop_waste.flow_mass_comp[j]
            )

        @self.Constraint(
            self.config.property_package.component_list, doc="Component mass balances"
        )
        def component_mass_balance(b, j):
            return (
                prop_in.flow_mass_comp[j]
                == prop_out.flow_mass_comp[j] + prop_waste.flow_mass_comp[j]
            )

        #
        @self.Constraint(doc="Outlet temperature equation")
        def isothermal_outlet(b):
            return prop_in.temperature == prop_out.temperature

        @self.Constraint(doc="Waste temperature equation")
        def isothermal_waste(b):
            return prop_in.temperature == prop_waste.temperature
        
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

        self.waste = Port(noruleinit=True, doc='Waste Port')
        self.waste.add(prop_waste.flow_vol, 'flow_vol')
        self.waste.add(prop_waste.conc_mass_comp, 'conc_mass')
        self.waste.add(prop_waste.temperature, 'temperature')
        self.waste.add(prop_waste.pressure, 'pressure')
        