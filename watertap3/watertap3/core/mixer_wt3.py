##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
##############################################################################
"""
General purpose mixer block for IDAES models
"""

import logging
import idaes.logger as idaeslog
from pyomo.environ import check_optimal_termination, Var, units as pyunits
from pyomo.common.config import ConfigBlock, ConfigValue, In
from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import (
    ConfigurationError,
)
from pyomo.network import Port
from idaes.core.solvers.get_solver import get_solver

__all__ = ["Mixer"]
solver = get_solver()

__author__ = "Kurban Sitterley"

# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("Mixer")
class MixerData(UnitModelBlockData):
    """ """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
**default** = False. Mixer1 blocks are always steady-state.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Mixer1 blocks do not contain holdup, thus this must be
False.""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for mixer",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
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
        "inlet_list",
        ConfigValue(
            domain=list,
            description="List of inlet names",
            doc="""A list containing names of inlets,
**default** - None.
**Valid values:** {
**None** - use num_inlets argument,
**list** - a list of names to use for inlets.}""",
        ),
    )
    CONFIG.declare(
        "num_inlets",
        ConfigValue(
            domain=int,
            description="Number of inlets to unit",
            doc="""Argument indicating number (int) of inlets to construct, not
used if inlet_list arg is provided,
**default** - None.
**Valid values:** {
**None** - use inlet_list arg instead, or default to 2 if neither argument
provided,
**int** - number of inlets to create (will be named with sequential integers
from 1 to num_inlets).}""",
        ),
    )
    CONFIG.declare(
        "construct_ports",
        ConfigValue(
            default=True,
            domain=In([True, False]),
            description="Construct inlet and outlet Port objects",
            doc="""Argument indicating whether model should construct Port objects
linked to all inlet states and the mixed state,
**default** - True.
**Valid values:** {
**True** - construct Ports for all states,
**False** - do not construct Ports.""",
        ),
    )

    def build(self):
        """ """
        # Call super.build()
        super().build()
        self.initialized = False
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True
        self.properties_out = prop_out = self.config.property_package.state_block_class(
            doc="Material properties of outlet stream", **tmp_dict
        )

        self.create_inlet_list()
        self.add_inlet_state_blocks()

        self.add_port_objects()

        # @self.Constraint(doc="Overall flow balance")
        # def flow_balance(b):
        #     return prop_out.flow_vol == sum(ib.flow_vol for ib in b.inlet_blocks)

        @self.Constraint(
            self.config.property_package.solute_set,
            doc="Component mass balances",
        )
        def component_mass_balance(b, j):
            return prop_out.flow_mass_comp[j] == sum(
                ib.flow_mass_comp[j] for ib in b.inlet_blocks
            )

        # @self.Constraint(doc="Outlet temperature equation")
        # def isothermal(b):
        #     return prop_out.temperature == sum(
        #         ib.temperature for ib in b.inlet_blocks
        #     ) / len(b.inlet_blocks)

        # @self.Constraint(doc="Outlet pressure equation")
        # def isobaric(b):
        #     return prop_out.pressure == sum(ib.pressure for ib in b.inlet_blocks) / len(
        #         b.inlet_blocks
        #     )

    def create_inlet_list(self):
        """
        Create list of inlet stream names based on config arguments.

        Returns:
            list of strings
        """
        if self.config.inlet_list is not None and self.config.num_inlets is not None:
            # If both arguments provided and not consistent, raise Exception
            if len(self.config.inlet_list) != self.config.num_inlets:
                raise ConfigurationError(
                    "{} Mixer1 provided with both inlet_list and "
                    "num_inlets arguments, which were not consistent ("
                    "length of inlet_list was not equal to num_inlets). "
                    "PLease check your arguments for consistency, and "
                    "note that it is only necessary to provide one of "
                    "these arguments.".format(self.name)
                )
        elif self.config.inlet_list is None and self.config.num_inlets is None:
            # If no arguments provided for inlets, default to num_inlets = 2
            self.config.num_inlets = 2

        # Create a list of names for inlet StateBlocks
        if self.config.inlet_list is not None:
            inlet_list = self.config.inlet_list
        else:
            inlet_list = [
                "inlet_" + str(n) for n in range(1, self.config.num_inlets + 1)
            ]
        self.inlet_list = inlet_list

    def add_inlet_state_blocks(self):
        """
        Construct StateBlocks for all inlet streams.

        Args:
            list of strings to use as StateBlock names

        Returns:
            list of StateBlocks
        """
        # Setup StateBlock argument dict
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True

        self.inlet_blocks = []

        for i in self.inlet_list:
            tmp_prop = self.config.property_package.state_block_class(
                doc=f"Material properties for {i}",
                **tmp_dict,
            )

            setattr(self, f"properties_{i}", tmp_prop)

            self.inlet_blocks.append(tmp_prop)

    def add_port_objects(self):
        """
        Adds Port objects if required.

        Args:
            a list of inlet StateBlock objects
            a mixed state StateBlock object

        Returns:
            None
        """
        if self.config.construct_ports is True:
            self.inlet_ports = list()

            for i, b in zip(self.inlet_list, self.inlet_blocks):
                tmp_port = Port(noruleinit=True, doc=f"{i.title()} Port")
                setattr(self, i, tmp_port)
                tmp_port.add(b.flow_mass_comp, "flow_mass")
                # tmp_port.add(b.pressure, "pressure")
                # tmp_port.add(b.temperature, "temperature")
                # tmp_port.add(b.conc_mass_comp, "conc_mass")
                # tmp_port.add(b.flow_vol, "flow_vol")
                self.inlet_ports.append(tmp_port)

            self.outlet = Port(noruleinit=True, doc="Outlet Port")
            self.outlet.add(self.properties_out.flow_mass_comp, "flow_mass")
            # self.outlet.add(self.properties_out.flow_vol, "flow_vol")
            # self.outlet.add(self.properties_out.conc_mass_comp, "conc_mass")
            # self.outlet.add(self.properties_out.temperature, "temperature")
            # self.outlet.add(self.properties_out.pressure, "pressure")

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        optarg={},
        solver=None,
        hold_state=True,
    ):
        """
        Initialization routine for mixer (default solver ipopt)

        Keyword Arguments:
            optarg : solver options dictionary object (default={})
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')
            hold_state : flag indicating whether the initialization routine
                     should unfix any state variables fixed during
                     initialization, **default** - False. **Valid values:**
                     **True** - states variables are not unfixed, and a dict of
                     returned containing flags for which states were fixed
                     during initialization, **False** - state variables are
                     unfixed after initialization by calling the release_state
                     method.

        Returns:
            If hold_states is True, returns a dict containing flags for which
            states were fixed during initialization.
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        if solver is None:
            opt = get_solver(solver, optarg)

        # Initialize inlet state blocks
        self.state_dict = state_dict = self.properties_out.define_port_members()
        self.flags = flags = {}
        self.flow_vol_list = list()
        self.state_args = state_args = dict(
            flow_vol=0,
            conc_mass_comp=dict(
                zip(
                    self.config.property_package.solute_set,
                    [0] * len(self.config.property_package.solute_set),
                )
            ),
        )
        for i, ib in zip(self.inlet_list, self.inlet_blocks):
            flags[i] = ib.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                hold_state=True,
                state_args=state_args,
            )
            init_log.info(f"Initialization for {i} on {self.name} Complete.")
            for k, v in state_dict.items():
                if k == "conc_mass_comp" and v.is_indexed():
                    for j in v.keys():
                        # self.conc_mass_dict[j].append(ib.conc_mass_comp[j].value)
                        self.state_args[k][j] += ib.conc_mass_comp[j].value
                elif k == "flow_vol":
                    self.state_args[k] += ib.flow_vol.value
                else:
                    pass

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=False,
            state_args=self.state_args,
        )
        init_log.info(f"Initialization for outlet on {self.name} Complete.")
        
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)

        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))
        for i, ib in zip(self.inlet_list, self.inlet_blocks):
            ib.release_state(flags[i], outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        self.initialized = True

    def calculate_scaling_facors(self):
        super().calculate_scaling_factors()

        # for b in self.inlet_blocks:



