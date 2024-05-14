# “Institute for the Design of Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE Framework) Copyright (c) 2019, by the software owners:
# The Regents of the University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University Research Corporation, et al.
# All rights reserved."
import os
import pandas as pd

from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.scaling import set_scaling_factor, get_scaling_factor

from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import Var, Param

__all__ = ["WT3UnitProcessBase"]
__author__ = "Kurban Sitterley"

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
chem_addition_params_file = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/data/chemical_addition_cost_curves.csv"
)

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

        self.initialized = False

    def handle_unit_params(self):
        for k, v in self.config.unit_params.items():
            if hasattr(self, k):
                p = getattr(self, k)
                if isinstance(p, Var):
                    p.fix(v)
                elif isinstance(p, Param):
                    p.set_value(v)
                else:
                    continue

    def get_chem_addition_params(self):
        """
        Get parameters used to cost chemical additions
        """

        if not self.chemical:
            raise ValueError("Must specify chemical name to get cost paramaters for.")

        df = pd.read_csv(chem_addition_params_file, index_col="chemical_name")
        self.chem_params_data = df.loc[self.chemical]

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        source_df = self.flowsheet().source_df

        if hasattr(self, "properties_in"):  # all
            if hasattr(self.properties_in, "flow_vol"):
                q_in = source_df.loc["flow"].value
                if get_scaling_factor(self.properties_in.flow_vol) is None:
                    set_scaling_factor(self.properties_in.flow_vol, 1 / q_in)
            if hasattr(self.properties_out, "flow_vol"):
                if get_scaling_factor(self.properties_out.flow_vol) is None:
                    set_scaling_factor(self.properties_out.flow_vol, 1 / q_in)
            for j, v in self.properties_in.flow_mass_comp.items():
                if get_scaling_factor(v) is None:
                    if j == "H2O":
                        flow_in = (
                            source_df.loc["flow"].value
                            * self.properties_in.dens_mass.value
                        )
                    else:
                        flow_in = source_df.loc[j].value * source_df.loc["flow"].value
                    set_scaling_factor(v, 1 / flow_in)
                if hasattr(self, "water_recovery"):
                    wrf = self.flowsheet().unit_water_recovery[self.unit_name]
                    if not isinstance(wrf, str) and wrf is not None:
                        pass
                    else:
                        wrf = 0.5
                    if get_scaling_factor(self.water_recovery) is None:
                        set_scaling_factor(self.water_recovery, 1 / wrf)
                    if (
                        get_scaling_factor(self.properties_out.flow_mass_comp[j])
                        is None
                    ):
                        if j == "H2O":
                            flow_out = wrf * flow_in
                        else:
                            try:
                                rf = self.flowsheet().unit_removal_fractions[
                                    self.unit_name
                                ][j]
                                flow_out = (1 - rf) * flow_in
                            except KeyError:
                                flow_out = flow_in
                        set_scaling_factor(
                            self.properties_out.flow_mass_comp[j], 1 / flow_out
                        )
                    if hasattr(self, "properties_waste"):
                        if (
                            get_scaling_factor(self.properties_waste.flow_mass_comp[j])
                            is None
                        ):
                            if j == "H2O":
                                flow_waste = (1 - wrf) * flow_in
                            else:
                                flow_waste = flow_in - flow_out
                            set_scaling_factor(
                                self.properties_waste.flow_mass_comp[j], 1 / flow_waste
                            )
                else:
                    flow_out = flow_in
                    if get_scaling_factor(
                        self.properties_out.flow_mass_comp[j] is None
                    ):
                        set_scaling_factor(
                            self.properties_out.flow_mass_comp[j], 1 / flow_out
                        )

                if hasattr(self.properties_in, "conc_mass_comp"):
                    if j == "H2O":
                        pass
                    else:
                        c_in = source_df.loc[j].value
                        try:
                            rf = self.flowsheet().unit_removal_fractions[
                                self.unit_name
                            ][j]
                            if rf == 0 or rf is None:
                                rf = 1e-5
                            c_out = (1 - rf) * c_in
                            c_waste = rf * c_in
                        except KeyError:
                            c_out = c_in
                            rf = 1
                        if (
                            get_scaling_factor(self.properties_in.conc_mass_comp[j])
                            is None
                        ):
                            set_scaling_factor(
                                self.properties_in.conc_mass_comp[j], 1 / c_in
                            )
                        if (
                            get_scaling_factor(self.properties_out.conc_mass_comp[j])
                            is None
                        ):
                            set_scaling_factor(
                                self.properties_out.conc_mass_comp[j], 1 / c_out
                            )
                        if hasattr(self, "removal_fraction"):
                            if get_scaling_factor(self.removal_fraction[j]) is None:
                                set_scaling_factor(self.removal_fraction[j], 1 / rf)
                            if hasattr(self, "properties_waste"):
                                if (
                                    get_scaling_factor(
                                        self.properties_waste.conc_mass_comp[j]
                                    )
                                    is None
                                ):
                                    set_scaling_factor(
                                        self.properties_waste.conc_mass_comp[j],
                                        1 / c_waste,
                                    )
