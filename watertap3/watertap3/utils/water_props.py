from idaes.core import (PhysicalParameterBlock, StateBlockData, declare_process_block_class)
# from idaes.core.components import Component
from idaes.core.base.components import Solute
from idaes.core.base.phases import LiquidPhase, AqueousPhase
from pyomo.environ import Set, units as pyunits
from pyomo.common.config import ConfigValue, In, Bool
# from . import generate_constituent_list
import pandas as pd
import numpy as np
import os

__all__ = ['WT3ParameterBlock',
           'WT3StateBlock']
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
crf_file = os.path.abspath(os.path.join(__location__, os.pardir)) + "/data/constituent_removal_factors.csv"

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
        '''
        Callable method for Block construction.
        '''
        super().build()

        self._state_block_class = WT3StateBlock

        self.Liq = AqueousPhase()
        self.component_list = Set(dimen=1)
        # train_constituent_list = self.generate_constituent_list()
        # self.generate_constituent_list()
        # self.train_constituent_list = ["tds", "toc"]

        for constituent_name in self.config.constituent_list:
            self.component_list.add(constituent_name)
            self.add_component(constituent_name, Solute())
            # setattr(self, constituent_name, Component())

    def generate_constituent_list(self):
        train = self.parent_block().train
        # getting the list of consituents with removal factors that are bigger than 0
        df = pd.read_csv(crf_file)
        df.case_study = np.where(df.case_study == 'default', train['case_study'], df.case_study)
        df = df[df.reference == train['reference']]
        df = df[df.case_study == train['case_study']]
        df = df[df.scenario == 'baseline']
        list1 = df[df.value >= 0].constituent.unique()
        list2 = self.parent_block().source_df.index.unique().to_list()

        self.parent_block().source_constituents = source_constituents = [x for x in list1 if x in list2]
        self.train_constituent_list = source_constituents

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units({
                'time': pyunits.s,
                'length': pyunits.m,
                'mass': pyunits.kg,
                'amount': pyunits.mol,
                'temperature': pyunits.K,
                # 'volume': pyunits.liter
                })

@declare_process_block_class("WT3StateBlock")
class WT3StateBlockData(StateBlockData):
    """
    This won't actually be used for most WaterTAP3 models, but is included to
    allow for future integration with ProteusLib and IDAES
    """

    def build(self):
        """
        Callable method for Block construction
        """
        return NotImplementedError()