# “Institute for the Design of Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE Framework) Copyright (c) 2019, by the software owners:
# The Regents of the University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University Research Corporation, et al.
# All rights reserved."
import idaes.logger as idaeslog
from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.util.config import is_physical_parameter_block
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import Var,Param, units as pyunits
from pyomo.network import Port

__all__ = ["WT3UnitProcessBase"]
__author__ = "Kurban Sitterley"

@declare_process_block_class("WT3UnitProcessBase")
class WT3UnitProcessBaseData(UnitModelBlockData):
    """
    Base unit model class for WaterTAP3.
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or 
        not,
        **default** = False. Equilibrium Reactors do not 
        support dynamic behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be 
        constructed or not.
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
            doc="""Property parameter object used to define property 
        calculations,
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
            description="Arguments to use for " "constructing property " "packages",
            doc="""A ConfigBlock with arguments to be 
        passed to a property block(s)
        and used when constructing these,
        **default** - None.
        **Valid values:** {
        see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "unit_params",
        ConfigValue(
            default=dict(),
            domain=dict,
            description="unit_params dict passed by the user",
            doc="""Dictionary of unit parameters that can be used to modify 
            values for unit model Vars from the input csv""",
        ),
    )

    def build(self):
        super().build()

    def handle_unit_params(self):
        for k, v in self.config.unit_params.items():
            if hasattr(self, k):
                p = getattr(self, k)
                if isinstance(p, Var):
                    p.fix(v)
                elif isinstance(p, Param):
                    p.set_value(v)
                else:
                    raise ValueError(f"{k} in unit_params for {self.name} is for a {type(p)}" 
                                     "but must be for a Var or Param."
                                     f"Remove {k} from unit_params on the input sheet.")
            