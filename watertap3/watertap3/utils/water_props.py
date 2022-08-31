from idaes.core import (PhysicalParameterBlock, StateBlockData, declare_process_block_class)
from idaes.core.components import Component
from idaes.core.phases import LiquidPhase
from pyomo.environ import units as pyunits

# from . import generate_constituent_list
import pandas as pd
import numpy as np

__all__ = ['WaterParameterBlock',
           'WaterStateBlock']


@declare_process_block_class("WaterParameterBlock")
class PhysicalParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class

    Define component and phase lists, along with base units
    """

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(PhysicalParameterData, self).build()

        self._state_block_class = WaterStateBlock

        self.Liq = LiquidPhase()
        # train_constituent_list = self.generate_constituent_list()
        self.generate_constituent_list()

        for constituent_name in self.train_constituent_list:
            setattr(self, constituent_name, Component())

    def generate_constituent_list(self):
        train = self.parent_block().train
        # getting the list of consituents with removal factors that are bigger than 0
        df = pd.read_csv('data/constituent_removal_factors.csv')
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
                'volume': pyunits.liter
                })

@declare_process_block_class("WaterStateBlock")
class WaterStateBlockData(StateBlockData):
    """
    This won't actually be used for most WaterTAP3 models, but is included to
    allow for future integration with ProteusLib and IDAES
    """

    def build(self):
        """
        Callable method for Block construction
        """
        return NotImplementedError()