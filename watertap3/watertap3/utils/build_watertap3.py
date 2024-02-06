import ast
import os

import numpy as np
import pandas as pd

from pyomo.environ import Block, ConcreteModel, Var
from pyomo.network import Arc
from idaes.core import FlowsheetBlock

from . import module_import
from .source_wt3 import Source
from watertap3.utils import Mixer, Splitter, SplitterBinary, financials
from .water_props import WaterParameterBlock
from watertap3.wt_units.wt_unit_pt import WT3UnitProcessPT
from watertap3.wt_units.wt_unit_siso import WT3UnitProcessSISO

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

treatment_train_setup_file = os.path.abspath(os.path.join(__location__, os.pardir)) + "/data/treatment_train_setup.csv"
case_study_water_sources_file = os.path.abspath(os.path.join(__location__, os.pardir)) + "/data/case_study_water_sources.csv"



__all__ = [
           'watertap3_setup',
           'get_case_study',
           'add_unit_process', 
           'add_water_source',
           'get_pfd_dict',
           'generate_constituent_list',
           'create_wt3_unit',
           'create_arcs',
           'create_arc_dict',
           'check_split_mixer_need',
           'create_mixers',
           'create_splitters',
           'add_waste_streams'
           ]

def watertap3_setup(dynamic=False, case_study=None, reference='nawi', scenario='baseline',
                   source_reference=None, source_case_study=None, source_scenario=None, 
                   new_df_units=None, print_it=True, ro_bounds='swro', desired_recovery=1):
    '''
    Initial setup of WaterTAP3 model. 

    Create flowsheet and read in basic information about model (water sources, units in treatment train)
    '''

    def get_source(reference, water_type, case_study, scenario):
        '''
        Read in water source data. 
        '''
        df = pd.read_csv(case_study_water_sources_file, index_col='variable')
        if 'test' in case_study:
            case_study = 'test'
        try:
            source_df = df[((df.case_study == case_study) & (df.water_type == water_type) & \
                (df.reference == reference) & (df.scenario == scenario))].copy()
            source_flow = source_df.loc['flow'].value
        except:
            source_df = df[((df.case_study == case_study) & (df.water_type == water_type) & \
                (df.reference == reference) & (df.scenario == 'baseline'))].copy()
            source_flow = source_df.loc['flow'].value
        source_df.drop(source_df[source_df.index == 'flow'].index, inplace=True)
        return source_flow, source_df
    
    case_study_print = case_study.replace('_', ' ').swapcase()
    scenario_print = scenario.replace('_', ' ').swapcase()

    m_name = f'{case_study_print}: {scenario_print}'
    m = ConcreteModel(name=m_name)

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.train = {
            'case_study': case_study,
            'reference': reference,
            'scenario': scenario
            }

    m.fs.case_study = case_study
    m.fs.scenario = scenario
    m.fs.reference = reference
    m.fs.ro_bounds = ro_bounds
    m.fs.desired_recovery = desired_recovery
    if source_reference is None:
        source_reference = reference
    if source_case_study is None:
        source_case_study = case_study
    if source_scenario is None:
        source_scenario = scenario

    df = pd.read_csv(treatment_train_setup_file) # Read in treatment train input sheet.

    water_type_list = []
    if new_df_units is not None:
        m.fs.df_units = new_df_units.copy()
    else:
        m.fs.df_units = df[((df.Reference == reference) & \
            (df.Scenario == scenario) & (df.CaseStudy == case_study))].copy()


    if 'binary_train' in m.fs.df_units.columns:
        cols = ['CaseStudy', 'binary_train', 'Reference', 'Scenario', 'Unit', \
            'Type', 'UnitName', 'ToUnitName', 'FromPort', 'Parameter']
    else:
        cols = ['CaseStudy', 'Reference', 'Scenario', 'Unit', 'Type', \
            'UnitName', 'ToUnitName', 'FromPort', 'Parameter']
    
    m.fs.df_units = m.fs.df_units[cols].copy()
    
    m.fs.has_ro = False
    m.fs.has_ix = False
    if 'ion_exchange' in m.fs.df_units.Unit:
        m.fs.has_ix = True
    if 'reverse_osmosis' in m.fs.df_units.Unit:
        m.fs.has_ro = True

    for i in m.fs.df_units[m.fs.df_units.Type == 'intake'].index:
        temp_dict = ast.literal_eval(m.fs.df_units[m.fs.df_units.Type == 'intake'].loc[i]['Parameter'])
        for water_type in temp_dict['water_type']:
            water_type_list.append(water_type)

    if len(water_type_list) == 1:
        water_type_list = water_type_list[0]


    m.fs.source_water = {
            'case_study': source_case_study,
            'reference': source_reference,
            'scenario': source_scenario,
            'water_type': water_type_list
            }

    flow_dict = {}

    if isinstance(m.fs.source_water['water_type'], list):
        m.fs.source_df = pd.DataFrame()
        for water_type in m.fs.source_water['water_type']:
            source_flow, source_df = get_source(m.fs.source_water['reference'],
                                                water_type,
                                                m.fs.source_water['case_study'],
                                                m.fs.source_water['scenario'])
            flow_dict[water_type] = source_flow
            m.fs.source_df = m.fs.source_df.append(source_df)


    else:
        source_flow, source_df = get_source(m.fs.source_water['reference'],
                                            m.fs.source_water['water_type'],
                                            m.fs.source_water['case_study'],
                                            m.fs.source_water['scenario'])
        flow_dict[m.fs.source_water['water_type']] = source_flow
        m.fs.source_df = source_df

    
    if print_it:
        print(f'\nCase Study = {case_study_print}'
                f'\nScenario = {scenario_print}')
                
    m.fs.flow_in_dict = flow_dict

    return m

def get_case_study(m, new_df_units=None, print_it=True):
    '''
    Function to add constituents and unit processes to flowsheet and connect 
    all ports.
    '''
    if new_df_units is not None:
        m.fs.df_units = new_df_units
        m.fs.pfd_dict = get_pfd_dict(new_df_units)
        m.fs.new_case_study = True
    else:
        m.fs.pfd_dict = get_pfd_dict(m.fs.df_units)
        m.fs.new_case_study = False

    pfd_dict = m.fs.pfd_dict
    financials.get_system_specs(m.fs)

    m.fs.water = WaterParameterBlock()

    if print_it:
        print('\n=========================ADDING UNIT PROCESSES=========================')
        for unit_process_name in pfd_dict.keys():
            unit = unit_process_name.replace('_', ' ').swapcase()
            unit_process_type = pfd_dict[unit_process_name]['Unit']
            unit_process_kind = pfd_dict[unit_process_name]['Type']
            print(f'{unit}')
            m = add_unit_process(m=m,
                    unit_process_name=unit_process_name,
                    unit_process_type=unit_process_type,
                    unit_process_kind=unit_process_kind)
        print('=======================================================================\n')
    else:
        for unit_process_name in pfd_dict.keys():
            unit = unit_process_name.replace('_', ' ').swapcase()
            unit_process_type = pfd_dict[unit_process_name]['Unit']
            unit_process_kind = pfd_dict[unit_process_name]['Type']
            m = add_unit_process(m=m,
                    unit_process_name=unit_process_name,
                    unit_process_type=unit_process_type,
                    unit_process_kind=unit_process_kind)

    # create a dictionary with all the arcs in the network based on the pfd_dict
    m, arc_dict, arc_i = create_arc_dict(m, pfd_dict, m.fs.flow_in_dict)
    m.fs.arc_dict = arc_dict

    # gets list of unit processes and ports that need either a splitter or mixer 
    splitter_list, mixer_list = check_split_mixer_need(arc_dict)
    m.fs.splitter_list = splitter_list
    m.fs.mixer_list = mixer_list
    # add the mixers if needed, and add the arcs around the mixers to the arc dictionary
    m, arc_dict, mixer_i, arc_i = create_mixers(m, mixer_list, arc_dict, arc_i)

    # add the splitters if needed, and add the arcs around the splitters to the arc dictionary
    m, arc_dict, splitter_i, arc_i = create_splitters(m, splitter_list, arc_dict, arc_i)
    m.fs.splitter_i = splitter_i
    # add the arcs to the model
    m = create_arcs(m, arc_dict)
    # add the waste arcs to the model
    m, arc_i, mixer_i = add_waste_streams(m, arc_i, pfd_dict, mixer_i)

    return m

def add_unit_process(m=None, unit_process_name=None,
        unit_process_type=None, unit_process_kind=None):

    up_module = module_import.get_module(unit_process_type)

    unit_params = m.fs.pfd_dict[unit_process_name]['Parameter']

    if unit_process_type == 'basic_unit':
        setattr(m.fs, unit_process_name,
            up_module.UnitProcess(default={'property_package': m.fs.water}))
        basic_unit_name = unit_params['unit_process_name']
        m = create_wt3_unit(m, basic_unit_name, unit_process_name)

    else:
        setattr(m.fs, unit_process_name,
            up_module.UnitProcess({"property_package": m.fs.water}))
        if not isinstance(getattr(m.fs, unit_process_name), WT3UnitProcessPT): 
            m = create_wt3_unit(m, unit_process_type, unit_process_name)

    unit = getattr(m.fs, unit_process_name)
    unit.chem_dict = {}
    unit.costing = Block()
    unit.tpec_tic = Var(
        bounds=(0, None),
        initialize=1,
        doc='Capital factor (TPEC or TIC)')
    unit.tpec_tic.fix(1)

    unit.unit_type = unit_process_type
    unit.unit_name = unit_process_name
    unit.unit_pretty_name = unit_process_name.replace('_', ' ').title().replace('Ro', 'RO').replace('Zld', 'ZLD').replace('Aop', 'AOP'). \
        replace('Uv', 'UV').replace('And', '&').replace('Sw', 'SW').replace('Gac', 'GAC').replace('Ph', 'pH').replace('Bc', 'BC'). \
        replace('Wwtp', 'WWTP').replace('Pac', 'PAC').replace('Co2', 'CO2').replace('Kmno4', 'KMnO4')
    unit.unit_kind = unit_process_kind
    if isinstance(unit_params, float):
        unit_params = {}
    unit.unit_params = unit_params
    unit.get_costing()

    return m

def add_water_source(m=None, source_name=None, water_type=None, flow=None, link_to=None):
    setattr(m.fs, source_name, Source(default={'property_package': m.fs.water}))
    getattr(m.fs, source_name).set_source()
    getattr(m.fs, source_name).flow_vol_in.fix(flow)
    temp_source_df = m.fs.source_df[m.fs.source_df.water_type == water_type].copy()
    train_constituent_list = list(getattr(m.fs, source_name).config.property_package.component_list)
    for constituent_name in train_constituent_list:
        if constituent_name in temp_source_df.index:
            conc = temp_source_df.loc[constituent_name].value
            getattr(m.fs, source_name).conc_mass_in[:, constituent_name].fix(conc)
        else:
            getattr(m.fs, source_name).conc_mass_in[:, constituent_name].fix(0)

    getattr(m.fs, source_name).pressure_in.fix(1)
    return m

def create_wt3_unit(m, unit_process_type, unit_process_name):

    def get_removal_factors():
        train = m.fs.train
        df = pd.read_csv('data/constituent_removal_factors.csv')
        const_df = df[((df.unit_process == unit_process_type) & (df.scenario == 'baseline') & \
            (df.reference == train['reference']))].copy()
        const_df = const_df[(const_df.case_study == train['case_study']) | \
            (const_df.case_study == 'default')].copy()
        constituent_list = getattr(m.fs, unit_process_name).config.property_package.component_list
        removal_dict = {}
        for constituent in constituent_list:
            if constituent not in const_df.constituent.unique():
                continue
            if const_df[((const_df.case_study == train['case_study']) & \
                (const_df.constituent == constituent))].empty:
                rf = const_df[((const_df.case_study == 'default') & \
                    (const_df.constituent == constituent))].value.iloc[0]
                removal_dict[constituent] = rf
            else:
                rf = const_df[((const_df.case_study == train['case_study']) & \
                    (const_df.constituent == constituent))].value.iloc[0]
                removal_dict[constituent] = rf

        return removal_dict

    wr_df = pd.read_csv('data/water_recovery_factors.csv')
    case_study_name = m.fs.train['case_study']
    scenario = m.fs.train['scenario']
    unit = getattr(m.fs, unit_process_name)


    cases = wr_df[wr_df.unit_process == unit_process_type].case_study.to_list()
    scenarios = wr_df[wr_df.unit_process == unit_process_type].scenario.to_list()
    default_df = wr_df[((wr_df.unit_process == unit_process_type) & \
        (wr_df.case_study == 'default'))].recovery
    tups = zip(cases, scenarios)

    if not isinstance(unit, WT3UnitProcessSISO):

        if (case_study_name, scenario) in tups:
            case_study_df = wr_df[((wr_df.unit_process == unit_process_type) & \
                (wr_df.case_study == case_study_name) & (wr_df.scenario == scenario))]
            if 'calculated' not in case_study_df.recovery.max():
                flow_recovery_factor = float(case_study_df.recovery)
                getattr(m.fs, unit_process_name).water_recovery.fix(flow_recovery_factor)
        else:
            if default_df.empty:
                raise TypeError(f'There is no default water recovery for {unit_process_type}.\n'
                                'Check that there is an entry for this unit in water_recovery.csv')
            if 'calculated' not in default_df.max():
                flow_recovery_factor = float(default_df)
                getattr(m.fs, unit_process_name).water_recovery.fix(flow_recovery_factor)

    train_constituent_removal_factors = \
        get_removal_factors()

    for constituent_name in getattr(m.fs, unit_process_name).config.property_package.component_list:
        
        if constituent_name in train_constituent_removal_factors.keys():
            unit.removal_fraction[:, constituent_name].fix(train_constituent_removal_factors[constituent_name])
        elif isinstance(unit, WT3UnitProcessSISO):
            unit.removal_fraction[:, constituent_name].fix(0)
        else:
            unit.removal_fraction[:, constituent_name].fix(1E-5)
    return m


def generate_constituent_list(m):
        train = m.fs.train
        # getting the list of consituents with removal factors that are bigger than 0
        df = pd.read_csv('data/constituent_removal_factors.csv')
        df.case_study = np.where(df.case_study == 'default', train['case_study'], df.case_study)
        df = df[df.reference == train['reference']]
        df = df[df.case_study == train['case_study']]
        df = df[df.scenario == 'baseline']
        list1 = df[df.value >= 0].constituent.unique()
        list2 = m.fs.source_df.index.unique().to_list()

        m.fs.source_constituents = source_constituents = [x for x in list1 if x in list2]

        return source_constituents

def get_pfd_dict(df_units):
    ### create pfd_dictionary for treatment train
    pfd_dict = df_units.set_index('UnitName').T.to_dict()
    for key in pfd_dict.keys():
        # parameter from string to dict
        if pfd_dict[key]['Parameter'] is not np.nan:
            pfd_dict[key]['Parameter'] = ast.literal_eval(pfd_dict[key]['Parameter'])
        if pfd_dict[key]['ToUnitName'] is not np.nan:
            if ',' in pfd_dict[key]['ToUnitName']:
                pfd_dict[key]['ToUnitName'] = pfd_dict[key]['ToUnitName'].split(',')
                pfd_dict[key]['FromPort'] = pfd_dict[key]['FromPort'].split(',')

    return pfd_dict

def create_arcs(m, arc_dict):
    for key in arc_dict.keys():
        source = arc_dict[key][0]
        source_port = arc_dict[key][1]
        outlet = arc_dict[key][2]
        outlet_port = arc_dict[key][3]
        setattr(m.fs, (f'arc{key}'), Arc(source=getattr(getattr(m.fs, source), source_port),
                                         destination=getattr(getattr(m.fs, outlet), outlet_port)))
    return m

def create_arc_dict(m, pfd_dict, flow):
    # create arc dictionary, add sources, add source to inlet arcs
    arc_dict = {}
    arc_i = 1
    for key in pfd_dict.keys():
        # if the unit is an intake process
        if pfd_dict[key]['Type'] == 'intake':
            num_sources = len(pfd_dict[key]['Parameter']['water_type'])
            num_unique_sources = len(np.unique(pfd_dict[key]['Parameter']['water_type']))

            ### check if multiple sources with same name for 1 intake
            if num_sources != num_unique_sources:
                print('error: multiple sources with same name for 1 intake')

            for water_type in pfd_dict[key]['Parameter']['water_type']:
                source_name = water_type
                water_type = water_type
                source_flow = flow[source_name]

                m = add_water_source(m=m, source_name=source_name, water_type=water_type, flow=source_flow)

                arc_dict[arc_i] = [source_name, 'outlet', key, 'inlet']
                arc_i += 1

                # create arcs *for single streams* from .csv table.
    for key in pfd_dict.keys():
        if pfd_dict[key]['FromPort'] is not np.nan:
            if isinstance(pfd_dict[key]['FromPort'], list):
                for port_i in range(0, len(pfd_dict[key]['FromPort'])):
                    arc_dict[arc_i] = [key, pfd_dict[key]['FromPort'][port_i],
                                       pfd_dict[key]['ToUnitName'][port_i], 'inlet']
                    arc_i += 1
            else:
                arc_dict[arc_i] = [key, pfd_dict[key]['FromPort'], pfd_dict[key]['ToUnitName'], 'inlet']
                arc_i += 1

    return m, arc_dict, arc_i

def check_split_mixer_need(arc_dict):
    # check if a mixer or splitter is needed
    mixer_list = []
    splitter_list = []
    unique_name_list1 = []
    unique_name_list2 = []

    for key in arc_dict.keys():
        # FOR SPLITTER
        if [arc_dict[key][0], arc_dict[key][1]] not in unique_name_list1:
            unique_name_list1.append([arc_dict[key][0], arc_dict[key][1]])
        else:
            if [arc_dict[key][0], arc_dict[key][1]] not in splitter_list:
                splitter_list.append([arc_dict[key][0], arc_dict[key][1]])

        # FOR MIXER    
        if [arc_dict[key][2], arc_dict[key][3]] not in unique_name_list2:
            unique_name_list2.append([arc_dict[key][2], arc_dict[key][3]])
        else:
            if [arc_dict[key][2], arc_dict[key][3]] not in mixer_list:
                mixer_list.append([arc_dict[key][2], arc_dict[key][3]])
    return splitter_list, mixer_list

def create_mixers(m, mixer_list, arc_dict, arc_i):
    mixer_i = 1
    inlet_i = 1
    for j in mixer_list:
        inlet_list = []
        mixer_name = f'mixer{mixer_i}'
        for key in list(arc_dict.keys()):
            if ((arc_dict[key][2] == j[0]) & (arc_dict[key][3] == j[1])):

                # inlet list for when mixer is added to model
                inlet_name = f'inlet{inlet_i}'
                inlet_list.append(inlet_name)
                inlet_i += 1

                # add new arc to arc dict
                arc_dict[arc_i] = [arc_dict[key][0], arc_dict[key][1], mixer_name, inlet_name]
                arc_i += 1

                # delete from arc dict
                del arc_dict[key]

        # add mixer to model with inlet list
        setattr(m.fs, mixer_name,
                Mixer(default={'property_package': m.fs.water, 'inlet_list': inlet_list}))

        # arc from mixer outlet to node
        arc_dict[arc_i] = [mixer_name, 'outlet', j[0], j[1]]
        arc_i += 1
        mixer_i += 1

    return m, arc_dict, mixer_i, arc_i

def create_splitters(m, splitter_list, arc_dict, arc_i):
    splitter_i = 1
    outlet_i = 1
    m.fs.unit_options = unit_options = {}
    m.fs.all_splitters = all_splitters = {}
    if not splitter_list:
        m.fs.choose = False
    for j in splitter_list:
        splitter_unit = j[0]
        splitter_port = j[1]
        outlet_i = 1
        outlet_list = []
        outlet_list_up = m.fs.outlet_list_up = {}
        splitter_name = f'splitter{splitter_i}'
        all_splitters[splitter_name] = {'from_unit': splitter_unit, 
                                        'to_units': list(zip(m.fs.pfd_dict[splitter_unit]['ToUnitName'], 
                                        m.fs.pfd_dict[splitter_unit]['FromPort'])), 
                                        # 'split_fraction': m.fs.pfd_dict[splitter_unit]['Parameter']['split_fraction'], 
                                        'indicator': False}
        for key in list(arc_dict.keys()):
            if ((arc_dict[key][0] == splitter_unit) & (arc_dict[key][1] == splitter_port)):
                split_dict = m.fs.split_dict = {}
                w = 0
                for uname in m.fs.pfd_dict[splitter_unit]['ToUnitName']:
                    if m.fs.pfd_dict[splitter_unit]['FromPort'][w] == 'outlet':
                        if 'split_fraction' in m.fs.pfd_dict[splitter_unit]['Parameter']:
                            all_splitters[splitter_name]['split_fraction']  = \
                                m.fs.pfd_dict[splitter_unit]['Parameter']['split_fraction']
                            split_dict[uname] = m.fs.pfd_dict[splitter_unit]['Parameter']['split_fraction'][w]
                            w += 1
                            split_fractions = m.fs.pfd_dict[splitter_unit]['Parameter']['split_fraction']
                            if all(split == 1 for split in split_fractions):
                                m.fs.choose = True
                                unit_options[splitter_i] = {splitter_unit: list(split_dict.keys())}
                            else:
                                m.fs.choose = False

                            # outlet list for when splitter is added to model
                outlet_name = f'outlet_{outlet_i}'
                outlet_list.append(outlet_name)
                outlet_i += 1

                # add new arc to arc dict
                arc_dict[arc_i] = [splitter_name, outlet_name, arc_dict[key][2], arc_dict[key][3]]
                arc_i += 1

                unit_hold = arc_dict[key][2]

                if arc_dict[key][2] not in split_dict.keys():
                    for l in m.fs.arc_dict.keys():
                        if arc_dict[key][2] == m.fs.arc_dict[l][0]:
                            unit_hold = m.fs.arc_dict[l][2]

                if len(split_dict) > 0:
                    outlet_list_up[outlet_name] = split_dict[unit_hold]
                else:
                    outlet_list_up[outlet_name] = "NA"
                # delete from arc dict
                del arc_dict[key]

        # add splitter to model with outlet list
        if m.fs.choose:
            all_splitters[splitter_name]['indicator'] = True
            setattr(m.fs, splitter_name, SplitterBinary(default={'property_package': m.fs.water}))
            this_splitter = getattr(m.fs, splitter_name)
            this_splitter.split_dict = split_dict
            this_splitter._split_from_unit = splitter_unit
            this_splitter.get_split(split_dict=split_dict)
            # setattr(this_splitter, f'unit_list', unit_options[splitter_i])
        else:
            setattr(m.fs, splitter_name, Splitter(default={'property_package': m.fs.water}))
            setattr(m.fs, f'{splitter_name}_outlet_list', outlet_list_up)
            unit_params = m.fs.pfd_dict[splitter_unit]['Parameter']
            getattr(m.fs, splitter_name).outlet_list = outlet_list_up
            # could just have self call the split list directly without reading in unit params. same for all 
            getattr(m.fs, splitter_name).get_split(outlet_list_up=outlet_list_up, unit_params=unit_params)

        # arc from mixer outlet to node
        arc_dict[arc_i] = [splitter_unit, splitter_port, splitter_name, 'inlet']
        arc_i += 1
        splitter_i += 1

    return m, arc_dict, splitter_i, arc_i

def add_waste_streams(m, arc_i, pfd_dict, mixer_i):
    # get number of units going to automatic waste disposal units
    i = 0
    waste_inlet_list = []
    unit_list = []
    for key in m.fs.pfd_dict.keys():
        unit_list.append(m.fs.pfd_dict[key]['Unit'])
        if 'surface_discharge' == m.fs.pfd_dict[key]['Unit']:
            sd_name = key
    if 'surface_discharge' in unit_list:
        for b_unit in m.fs.component_objects(Block, descend_into=False):
            if hasattr(b_unit, 'waste'):
                if len(getattr(b_unit, 'waste').arcs()) == 0:
                    if str(b_unit)[3:] in list(pfd_dict.keys()):  #
                        if pfd_dict[str(b_unit)[3:]]['Type'] == 'treatment':
                            i += 1
                            waste_inlet_list.append((f'inlet{i}'))

        if len(waste_inlet_list) > 1:
            i = 0
            waste_mixer = f'mixer{mixer_i}'
            m.fs.water_mixer_name = waste_mixer  # used for displaying train. not used for model
            setattr(m.fs, waste_mixer,
                    Mixer(default={'property_package': m.fs.water, 'inlet_list': waste_inlet_list}))
            for b_unit in m.fs.component_objects(Block, descend_into=False):
                if hasattr(b_unit, 'waste'):
                    if len(getattr(b_unit, 'waste').arcs()) == 0:
                        if str(b_unit)[3:] in list(pfd_dict.keys()):
                            if pfd_dict[str(b_unit)[3:]]['Type'] == 'treatment':
                                i += 1
                                setattr(m.fs, (f'arc{arc_i}'), \
                                    Arc(source=getattr(b_unit, 'waste'),
                                    destination=getattr(getattr(m.fs, waste_mixer), f'inlet{i}')))
                                arc_i += 1

            # add connection for waste mixer to surface discharge -->
            if 'surface_discharge' in unit_list:
                setattr(m.fs, (f'arc{arc_i}'), \
                    Arc(source=getattr(m.fs, waste_mixer).outlet,
                    destination=getattr(m.fs, sd_name).inlet))
                arc_i += 1
        return m, arc_i, mixer_i

    else:
        return m, arc_i, mixer_i