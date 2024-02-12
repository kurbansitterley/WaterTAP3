from idaes.core import (
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    declare_process_block_class,
)

# from idaes.core.components import Component
from idaes.core.base.components import Solute, Solvent
from idaes.core.base.phases import LiquidPhase, AqueousPhase
from pyomo.environ import Set, units as pyunits
from pyomo.common.config import ConfigValue, In, Bool
from idaes.core.util.constants import Constants

# from . import generate_constituent_list
import pandas as pd
import numpy as np
import os
import idaes.logger as idaeslog
from idaes.core.util.misc import add_object_reference
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_unfixed_variables,
)
from pyomo.environ import (
    Param,
    PositiveReals,
    units as pyunits,
    Var,
    Constraint,
    Suffix,
    value,
    check_optimal_termination,
)

__all__ = ["WT3ParameterBlock", "WT3StateBlock"]
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
crf_file = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/data/constituent_removal_factors.csv"
)


@declare_process_block_class("WT3ParameterBlock")
class WT3ParameterBlockData(PhysicalParameterBlock):
    """
    Property Parameter Block Class

    Define component and phase lists, along with base units
    """

    CONFIG = PhysicalParameterBlock.CONFIG()
    CONFIG.declare(
        "constituent_list",
        ConfigValue(
            domain=list,
            description="Required argument. List of strings that specify names of constituents.",
        ),
    )

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        # self.config.constituent_list += ["H2O"]

        self._state_block_class = WT3StateBlock

        self.Liq = LiquidPhase()
        self.H2O = Solvent()
        # self.component_list = Set(dimen=1)
        # train_constituent_list = self.generate_constituent_list()
        # self.generate_constituent_list()
        # self.train_constituent_list = ["tds", "toc"]

        for j in self.config.constituent_list:
            # self.component_list.add(constituent_name)
            self.add_component(str(j), Solute())
            # setattr(self, constituent_name, Component())

        self.dens_mass = Param(
            initialize=1000,
            units=pyunits.kg / pyunits.m**3,
            mutable=True,
            doc="Mass density of pure water",
        )

        self.dens_mass_sw = Param(
            initialize=1030,
            units=pyunits.kg / pyunits.m**3,
            mutable=True,
            doc="Mass density of seawater",
        )

        self.visc_d = Param(
            initialize=0.001,
            units=pyunits.kg / pyunits.m / pyunits.s,
            mutable=True,
            doc="Dynamic viscosity of solution",
        )

        self.mw_sw = Param(
            initialize=58.4e-3,  # molecular weight of NaCl
            units=pyunits.kg / pyunits.mol,
            mutable=True,
            doc="Molecular weight of seawater",
        )

        # self.set_default_scaling("temperature", 1e-3)
        # self.set_default_scaling("pressure", 1e-5)
        self.set_default_scaling("pressure_osmotic", 1e-5)
        self.set_default_scaling("dens_mass", 1e-3)
        self.set_default_scaling("visc_d", 1e3)

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )
        obj.add_properties(
            {
                "flow_vol": {"method": None},
                "conc_mass_comp": {"method": None},
                "flow_mass_comp": {"method": "_flow_mass_comp"},
                "mass_frac_comp": {"method": "_mass_frac_comp"},
                "temperature": {"method": "_temperature"},
                "pressure": {"method": "_pressure"},
                "pressure_osmotic": {"method": "_pressure_osmotic"},
                "dens_mass": {"method": "_dens_mass"},
                "visc_d": {"method": "_visc_d"},
            }
        )


class _WT3StateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(
        self,
        state_args=None,
        state_vars_fixed=False,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Initialization routine for property package.

        Keyword Arguments:
        state_args : Dictionary with initial guesses for the state vars
                     chosen. Note that if this method is triggered
                     through the control volume, and if initial guesses
                     were not provied at the unit model level, the
                     control volume passes the inlet values as initial
                     guess.
        outlvl : sets output level of initialization routine
        state_vars_fixed: Flag to denote if state vars have already been
                          fixed.
                          - True - states have already been fixed and
                                   initialization does not need to worry
                                   about fixing and unfixing variables.
                         - False - states have not been fixed. The state
                                   block will deal with fixing/unfixing.
        optarg : solver options dictionary object (default=None, use
                 default solver options)
        solver : str indicating which solver to use during
                 initialization (default = None, use default solver)
        hold_state : flag indicating whether the initialization routine
                     should unfix any state variables fixed during
                     initialization (default=False).
                     - True - states varaibles are not unfixed, and
                             a dict of returned containing flags for
                             which states were fixed during
                             initialization.
                    - False - state variables are unfixed after
                             initialization by calling the
                             relase_state method

        Returns:
            If hold_states is True, returns a dict containing flags for
            which states were fixed during initialization.
        """

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")

        # Set solver and options
        opt = get_solver(solver, optarg)

        # Fix state variables
        flags = fix_state_vars(self, state_args)

        # initialize vars calculated from state vars
        for k in self.keys():
            for j in self[k].params.component_list:
                if self[k].is_property_constructed("flow_mass_comp"):
                    if j == "H2O":
                        self[k].flow_mass_comp[j].set_value(
                            self[k].flow_vol * self[k].dens_mass
                        )
                    else:
                        self[k].flow_mass_comp[j].set_value(
                            self[k].flow_vol * self[k].conc_mass_comp[j]
                        )

        # Check when the state vars are fixed already result in dof 0
        for k in self.keys():
            dof = degrees_of_freedom(self[k])
            if dof != 0:
                raise InitializationError(
                    "\nWhile initializing {sb_name}, the degrees of freedom "
                    "are {dof}, when zero is required. \nInitialization assumes "
                    "that the state variables should be fixed and that no other "
                    "variables are fixed. \nIf other properties have a "
                    "predetermined value, use the calculate_state method "
                    "before using initialize to determine the values for "
                    "the state variables and avoid fixing the property variables."
                    "".format(sb_name=self.name, dof=dof)
                )

        # # ---------------------------------------------------------------------
        skip_solve = True  # skip solve if only state variables are present
        for k in self.keys():
            if number_unfixed_variables(self[k]) != 0:
                skip_solve = False

        # if not skip_solve:
        #     # Initialize properties
        #     with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
        #         results = solve_indexed_blocks(opt, [self], tee=slc.tee)
        #         if not check_optimal_termination(results):
        #             raise InitializationError(
        #                 "The property package failed to solve during initialization."
        #             )
        #     init_log.info_high(
        #         "Property initialization: {}.".format(idaeslog.condition(results))
        #     )

        # ---------------------------------------------------------------------
        # If input block, return flags, else release state
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                self.release_state(flags)

    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        """
        Method to release state variables fixed during initialization.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")

        if flags is None:
            return

        # Unfix state variables
        revert_state_vars(self, flags)
        init_log.info("State Released.")


@declare_process_block_class("WT3StateBlock", block_class=_WT3StateBlock)
class WT3StateBlockData(StateBlockData):
    """
    This won't actually be used for most WaterTAP3 models, but is included to
    allow for future integration with ProteusLib and IDAES
    """

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)
        # self.H2O = Solvent()

        self.flow_vol = Var(
            initialize=1,
            domain=PositiveReals,
            units=pyunits.m**3 / pyunits.s,
            doc="Volumetric flow rate",
        )

        self.conc_mass_comp = Var(
            self.params.solute_set,
            initialize=1,
            domain=PositiveReals,
            units=pyunits.kg / pyunits.m**3,
            doc="Mass concentration of each solute",
        )

    # Other properties
    def _flow_mass_comp(self):

        self.flow_mass_comp = Var(
            self.params.component_list,
            initialize=1e2,
            domain=PositiveReals,
            doc="Mass flowrate of each component",
            units=pyunits.kg / pyunits.s,
        )

        def rule_flow_mass_comp(b, j):
            if j == "H2O":
                return b.flow_mass_comp[j] == b.flow_vol * b.dens_mass
            else:
                return b.flow_mass_comp[j] == b.flow_vol * b.conc_mass_comp[j]

        self.eq_flow_mass_comp = Constraint(
            self.params.component_list, rule=rule_flow_mass_comp
        )

    def _temperature(self):
        self.temperature = Var(
            initialize=298.15,
            bounds=(273.15, 373.15),
            units=pyunits.K,
            doc="Temperature",
        )

    def _pressure(self):
        self.pressure = Var(
            initialize=101325,
            bounds=(1e4, 8e6),
            units=pyunits.Pa,
            doc="Pressure",
        )

    def _mass_frac_comp(self):

        self.mass_frac_comp = Var(
            self.params.component_list,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Mass fraction",
        )

        def rule_mass_frac_comp(b, j):
            return b.mass_frac_comp[j] == b.flow_mass_comp[j] / sum(
                b.flow_mass_comp[j] for j in self.params.component_list
            )

        self.eq_mass_frac_comp = Constraint(
            self.params.component_list, rule=rule_mass_frac_comp
        )

    def _pressure_osmotic(self):
        self.pressure_osmotic = Var(
            initialize=1e6,
            bounds=(5e2, 5e7),
            units=pyunits.Pa,
            doc="Osmotic pressure",
        )

        self.osmotic_coefficient = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Osmotic coefficient",
        )

        self.osmotic_coeff_eq_A = Param(
            initialize=4.92,
            units=pyunits.dimensionless,
            doc="Osmotic coefficient relationship - A parameter",
        )
        self.osmotic_coeff_eq_B = Param(
            initialize=0.0889,
            units=pyunits.dimensionless,
            doc="Osmotic coefficient relationship - B parameter",
        )
        self.osmotic_coeff_eq_intercept = Param(
            initialize=0.918,
            units=pyunits.dimensionless,
            doc="Osmotic coefficient relationship - intercept",
        )

        def rule_osmotic_coefficient(b):
            return (
                b.osmotic_coefficient
                == b.osmotic_coeff_eq_A * b.mass_frac_comp["tds"] ** 2
                + b.osmotic_coeff_eq_B * b.mass_frac_comp["tds"]
                + b.osmotic_coeff_eq_intercept
            )

        self.eq_osmotic_coefficient = Constraint(rule=rule_osmotic_coefficient)

        def rule_osmotic_pressure(b):
            num_ions = 2
            molality_tds = (
                b.mass_frac_comp["tds"] / (1 - b.mass_frac_comp["tds"])
            ) / b.params.mw_sw
            temp = 293*pyunits.degK
            return b.pressure_osmotic == pyunits.convert(
                num_ions
                * b.osmotic_coefficient
                * molality_tds
                * b.params.dens_mass
                * Constants.gas_constant
                * temp,
                to_units=pyunits.Pa,
            )

        self.eq_osmotic_pressure = Constraint(rule=rule_osmotic_pressure)

    def _dens_mass(self):
        add_object_reference(self, "dens_mass", self.params.dens_mass)

    def _visc_d(self):
        add_object_reference(self, "visc_d", self.params.visc_d)

    def get_material_flow_terms(self, j):
        return self.flow_mass_comp[j]

    def get_enthalpy_flow_terms(self, p):
        raise NotImplementedError

    def get_material_density_terms(self, j):
        if j == "H2O":
            return self.dens_mass
        else:
            return self.conc_mass_comp[j]

    def get_energy_density_terms(self, p):
        raise NotImplementedError

    # def default_material_balance_type(self):
    #     return MaterialBalanceType.componentTotal

    # def default_energy_balance_type(self):
    #     return EnergyBalanceType.none

    def define_state_vars(self):
        return {"flow_vol": self.flow_vol, "conc_mass_comp": self.conc_mass_comp}

    def define_display_vars(self):
        return {
            "Volumetric Flowrate": self.flow_vol,
            "Mass Concentration": self.conc_mass_comp,
            # "Temperature": self.temperature,
        }

    # def get_material_flow_basis(self):
    #     return MaterialFlowBasis.mass

    def calculate_scaling_factors(self):

        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.flow_vol) is None:
            sf_Q = iscale.get_scaling_factor(self.flow_vol, default=1, warning=True)
            iscale.set_scaling_factor(self.flow_vol, sf_Q)

        for j, v in self.conc_mass_comp.items():
            sf_c = iscale.get_scaling_factor(self.conc_mass_comp[j])
            if sf_c is None:
                try:
                    sf_c = self.params.default_scaling_factor[("conc_mass_comp", j)]
                except KeyError:
                    iscale.set_scaling_factor(self.conc_mass_comp[j], 10)

        if self.is_property_constructed("flow_mass_comp"):
            for j, v in self.flow_mass_comp.items():
                if iscale.get_scaling_factor(v) is None:
                    if j == "H2O":
                        sf = value(self.flow_vol * self.dens_mass) ** -1
                    else:
                        sf = value(self.flow_vol * self.conc_mass_comp[j]) ** -1
                    iscale.set_scaling_factor(v, sf)
