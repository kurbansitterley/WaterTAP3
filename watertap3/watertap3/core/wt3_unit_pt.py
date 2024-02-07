
# “Institute for the Design of Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE Framework) Copyright (c) 2019, by the software owners:
# The Regents of the University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University Research Corporation, et al.
# All rights reserved."
import idaes.logger as idaeslog
from idaes.core import (UnitModelBlockData, declare_process_block_class, useDefault)
from idaes.core.util.config import is_physical_parameter_block
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import Var, units as pyunits
from pyomo.network import Port

__all__ = ['WT3UnitProcessPT']


@declare_process_block_class('WT3UnitProcessPT')
class WT3UnitProcessPTData(UnitModelBlockData):
    '''
    This class describes the rules for a zeroth-order model for a unit

    The Config Block is used to process arguments from when the model is
    instantiated. In IDAES, this serves two purposes:
         1. Allows us to separate physical properties from unit models
         2. Lets us give users options for configuring complex units
    The dynamic and has_holdup options are expected arguments which must exist
    The property package arguments let us define different sets of contaminants
    without needing to write a new model.
    '''

    CONFIG = ConfigBlock()
    CONFIG.declare('dynamic', ConfigValue(domain=In([False]), 
        default=False,
        description='Dynamic model flag - must be False',
        doc='''Indicates whether this model will be dynamic or 
        not,
        **default** = False. Equilibrium Reactors do not 
        support dynamic behavior.'''))
    CONFIG.declare('has_holdup', ConfigValue(default=False, 
        domain=In([False]),
        description='Holdup construction flag - must be False',
        doc='''Indicates whether holdup terms should be 
        constructed or not.
        **default** - False. Equilibrium reactors do not have defined volume, thus
        this must be False.'''))
    CONFIG.declare('property_package', ConfigValue(default=useDefault,
        domain=is_physical_parameter_block,
        description='Property package to use for control volume',
        doc='''Property parameter object used to define property 
        calculations,
        **default** - useDefault.
        **Valid values:** {
        **useDefault** - use default package from parent model or flowsheet,
        **PhysicalParameterObject** - a PhysicalParameterBlock object.}'''))
    CONFIG.declare('property_package_args', ConfigBlock(implicit=True,
        description='Arguments to use for '
                    'constructing property '
                    'packages',
        doc='''A ConfigBlock with arguments to be 
        passed to a property block(s)
        and used when constructing these,
        **default** - None.
        **Valid values:** {
        see property package for documentation.}'''))

    def build(self):
        super( ).build()
        # units_meta = self.config.property_package.get_metadata().get_derived_units
        time = self.flowsheet().config.time
        

        ## INLET
        self.flow_vol_in = Var(time,
            initialize=1,
            units=pyunits.m**3 / pyunits.s,
            bounds=(0, None),
            doc='Volumetric flowrate of water into unit')
        self.conc_mass_in = Var(time,
            self.config.property_package.component_list,
            initialize=1E-5,
            units=pyunits.kg / pyunits.m**3,
            doc='Mass concentration of species at inlet')
        self.temperature_in = Var(time,
            initialize=300,
            units=pyunits.degK,
            doc='Temperature at inlet')
        self.pressure_in = Var(time,
            initialize=1,
            units=pyunits.Pa,
            doc='Pressure at inlet')

        ## OUTLET
        self.flow_vol_out = Var(time,
            initialize=1,
            units=pyunits.m**3 / pyunits.s,
            doc='Volumetric flowrate of water out of unit')
        self.conc_mass_out = Var(time,
            self.config.property_package.component_list,
            initialize=0,
            units=pyunits.kg / pyunits.m**3,
            doc='Mass concentration of species at outlet')
        self.temperature_out = Var(time,
            initialize=300,
            units=pyunits.degK,
            doc='Temperature at outlet')
        self.pressure_out = Var(time,
            initialize=1,
            units=pyunits.Pa,
            doc='Pressure at outlet')
        self.deltaP_outlet = Var(time,
            initialize=1E-6,
            units=pyunits.Pa,
            doc='Pressure change between inlet and outlet')

        self.deltaP_outlet.fix(0)

        @self.Constraint(time, doc='Outlet pressure equation')
        def outlet_pressure_constraint(b, t):
            return (b.pressure_in[t] + b.deltaP_outlet[t] ==
                    b.pressure_out[t])

        @self.Constraint(time, doc='Overall flow balance')
        def flow_balance(b, t):
            return b.flow_vol_in[t] == b.flow_vol_out[t]

        @self.Constraint(time,
                         self.config.property_package.component_list,
                         doc='Component removal equation')
        def component_removal_equation(b, t, j):
            return (b.flow_vol_in[t] * b.conc_mass_in[t, j]  ==
                    b.flow_vol_out[t] * b.conc_mass_out[t, j])

        @self.Constraint(time, doc='Outlet temperature equation')
        def outlet_temperature_constraint(b, t):
            return b.temperature_in[t] == b.temperature_out[t]

        self.inlet = Port(noruleinit=True, doc='Inlet Port')
        self.inlet.add(self.flow_vol_in, 'flow_vol')
        self.inlet.add(self.conc_mass_in, 'conc_mass')
        self.inlet.add(self.temperature_in, 'temperature')
        self.inlet.add(self.pressure_in, 'pressure')

        self.outlet = Port(noruleinit=True, doc='Outlet Port')
        self.outlet.add(self.flow_vol_out, 'flow_vol')
        self.outlet.add(self.conc_mass_out, 'conc_mass')
        self.outlet.add(self.temperature_out, 'temperature')
        self.outlet.add(self.pressure_out, 'pressure')

        