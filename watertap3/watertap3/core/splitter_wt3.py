##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
##############################################################################


from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.util.config import is_physical_parameter_block
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import Constraint, NonNegativeReals, Var, units as pyunits
from pyomo.network import Port

module_name = "splitter_wt3"

__all__ = ["Splitter"]


@declare_process_block_class("Splitter")
class SplitterProcessData(UnitModelBlockData):
    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
        **default** = False. Equilibrium Reactors do not support dynamic behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
        **default** - False. Equilibrium reactors do not have defined volume, thus
        this must be False.""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
        **default** - useDefault.
        **Valid values:** {
        **useDefault** - use default package from parent model or flowsheet,
        **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
        and used when constructing these,
        **default** - None.
        **Valid values:** {
        see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "outlet_dict",
        ConfigValue(
            default=dict(),
            domain=dict,
            description="Dict of outlet names and corresponding split fractions",
        ),
    )

    def build(self):
        super().build()

        self.split_fraction_vars = []
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True
        self.properties_in = prop_in = self.config.property_package.state_block_class(
            doc="Material properties of inlet stream", **tmp_dict
        )

        self.inlet = Port(noruleinit=True, doc="Inlet Port")
        self.inlet.add(prop_in.flow_vol, "flow_vol")
        self.inlet.add(prop_in.conc_mass_comp, "conc_mass")
        self.inlet.add(prop_in.temperature, "temperature")
        self.inlet.add(prop_in.pressure, "pressure")

        for outlet, split in self.config.outlet_dict.items():
            tmp_port = Port(noruleinit=True, doc=f"{outlet.title()} Port")
            tmp_prop = self.config.property_package.state_block_class(
                doc=f"Material properties of {outlet} stream", **tmp_dict
            )
            tmp_split_var_name = f"split_fraction_{outlet}"
            tmp_split_var = Var(
                initialize=0.5,
                bounds=(0.001, 0.999),
                units=pyunits.dimensionless,
                doc=f"Split fraction for {outlet}",
            )
            tmp_split_var.construct()

            if split != "NA":
                tmp_split_var.fix(split)

            setattr(self, outlet, tmp_port)
            setattr(self, f"properties_{outlet}", tmp_prop)

            setattr(
                tmp_prop,
                tmp_split_var_name,
                tmp_split_var,
            )

            self.split_fraction_vars.append(tmp_split_var)

            tmp_port.add(tmp_prop.conc_mass_comp, "conc_mass")
            tmp_port.add(tmp_prop.flow_vol, "flow_vol")
            tmp_port.add(tmp_prop.temperature, "temperature")
            tmp_port.add(tmp_prop.pressure, "pressure")

            @tmp_prop.Constraint(
                self.config.property_package.component_list,
                doc=f"Component balance for {outlet}",
            )
            def component_mass_balance(b, j):
                return prop_in.conc_mass_comp[j] == b.conc_mass_comp[j]

            @tmp_prop.Constraint(doc=f"Flow split for {outlet}")
            def flow_split(b):
                return tmp_split_var * prop_in.flow_vol == b.flow_vol

            @tmp_prop.Constraint(doc=f"Isothermal for {outlet}")
            def isothermal_split(b):
                return b.temperature == prop_in.temperature

            @tmp_prop.Constraint(doc=f"Isobaric for {outlet}")
            def isobaric_split(b):
                return b.pressure == prop_in.pressure

        self.split_fraction_tot_ub = Constraint(
            expr=sum(self.split_fraction_vars) <= 1.025
        )
        self.split_fraction_tot_lb = Constraint(
            expr=sum(self.split_fraction_vars) >= 0.975
        )
