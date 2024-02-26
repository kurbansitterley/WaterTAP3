##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
##############################################################################
import idaes.logger as idaeslog
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.scaling import set_scaling_factor, get_scaling_factor
from pyomo.environ import Expression, check_optimal_termination, value, units as pyunits
from pyomo.network import Port
from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.util.config import is_physical_parameter_block
from pyomo.common.config import ConfigBlock, ConfigValue, In

module_name = "source_wt3"


@declare_process_block_class("Source")
class SourceData(UnitModelBlockData):
    """
    This class describes the rules for a zeroth-order model for a source.
    """

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

    def build(self):
        super().build()
        self.initialized = False
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False
        self.properties = prop = self.config.property_package.state_block_class(
            doc="Material properties of source", **tmp_dict
        )

        # prop.pressure.fix()
        prop.flow_vol
        self.flow_mgd = Expression(
            expr=pyunits.convert(
                self.properties.flow_vol, to_units=pyunits.Mgallons / pyunits.day
            )
        )

        self.outlet = Port(noruleinit=True, doc="Source Port")
        self.outlet.add(prop.flow_mass_comp, "flow_mass")
        # self.outlet.add(prop.flow_vol, 'flow_vol')
        # self.outlet.add(prop.conc_mass_comp, 'conc_mass')
        # self.outlet.add(prop.temperature, 'temperature')
        # self.outlet.add(prop.pressure, 'pressure')

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for initialization routines

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")

        if state_args is None:
            state_args = {}
            state_dict = self.properties.define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        self.properties.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )

        init_log.info("Initialization Complete")

        self.initialized = True

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
