##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
##############################################################################

from pyomo.network import Port
from pyomo.environ import check_optimal_termination, Var, Constraint, units as pyunits
from pyomo.common.config import ConfigBlock, ConfigValue, In

import idaes.logger as idaeslog
from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.solvers.get_solver import get_solver

module_name = "splitter_wt3"

__all__ = ["Splitter"]

__author__ = "Kurban Sitterley"


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

        self.initialized = False

        self.split_fraction_vars = []
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True
        self.properties_in = prop_in = self.config.property_package.state_block_class(
            doc="Material properties of inlet stream", **tmp_dict
        )

        self.inlet = Port(noruleinit=True, doc="Inlet Port")
        self.inlet.add(prop_in.flow_mass_comp, "flow_mass")
        # self.inlet.add(prop_in.flow_vol, "flow_vol")
        # self.inlet.add(prop_in.conc_mass_comp, "conc_mass")
        # self.inlet.add(prop_in.temperature, "temperature")
        # self.inlet.add(prop_in.pressure, "pressure")
        self.outlet_blocks = list()
        self.outlet_ports = list()
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
            self.outlet_blocks.append(tmp_prop)
            self.outlet_ports.append(tmp_port)

            tmp_port.add(tmp_prop.flow_mass_comp, "flow_mass")
            # tmp_port.add(tmp_prop.conc_mass_comp, "conc_mass")
            # tmp_port.add(tmp_prop.flow_vol, "flow_vol")
            # tmp_port.add(tmp_prop.temperature, "temperature")
            # tmp_port.add(tmp_prop.pressure, "pressure")

            @tmp_prop.Constraint(
                self.config.property_package.component_list,
                doc=f"Component balance for {outlet}",
            )
            def component_mass_balance(b, j):
                # return prop_in.conc_mass_comp[j] == b.conc_mass_comp[j]
                return tmp_split_var * prop_in.flow_mass_comp[j] == b.flow_mass_comp[j]

            # @tmp_prop.Constraint(doc=f"Flow split for {outlet}")
            # def flow_split(b):
            #     return tmp_split_var * prop_in.flow_vol == b.flow_vol

            # @tmp_prop.Constraint(doc=f"Isothermal for {outlet}")
            # def isothermal_split(b):
            #     return b.temperature == prop_in.temperature

            # @tmp_prop.Constraint(doc=f"Isobaric for {outlet}")
            # def isobaric_split(b):
            #     return b.pressure == prop_in.pressure

        self.split_fraction_tot_ub = Constraint(
            expr=sum(self.split_fraction_vars) <= 1.025
        )
        self.split_fraction_tot_lb = Constraint(
            expr=sum(self.split_fraction_vars) >= 0.975
        )

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
        self.state_dict = state_dict = self.properties_in.define_port_members()

        for ob, (outlet, split) in zip(
            self.outlet_blocks, self.config.outlet_dict.items()
        ):
            tmp_state_args = dict(
                flow_vol=0,
                conc_mass_comp=dict(
                    zip(
                        self.config.property_package.solute_set,
                        [0] * len(self.config.property_package.solute_set),
                    )
                ),
            )
            for k, v in state_dict.items():
                if k == "conc_mass_comp" and v.is_indexed():
                    for j in v.keys():
                        tmp_state_args[k][j] = self.properties_in.conc_mass_comp[
                            j
                        ].value
                elif k == "flow_vol":
                    tmp_state_args[k] = v.value * split
                else:
                    pass
                    # raise InitializationError?
            ob.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                hold_state=False,
                state_args=tmp_state_args,
            )
            init_log.info(f"Initialization for {outlet} on {self.name} complete.")

        flags = self.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
            state_args=state_args,
        )
        init_log.info(f"Initialization for inlet on {self.name} complete.")
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)

        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))
        self.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        self.initialized = True
