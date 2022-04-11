import ast
import logging
from re import template
import warnings
import pyomo.util.infeasible as infeas
import numpy as np
import pandas as pd
from idaes.core.util.model_statistics import degrees_of_freedom
import os
from pyomo.environ import Constraint, Objective, SolverFactory, \
    TransformationFactory, value, units as pyunits
from pyomo.network import SequentialDecomposition, Arc
from pyomo.network.port import SimplePort

from . import financials
from .build_watertap3 import *
from .post_processing import get_results_table
from pyomo.gdp import *

from watertap3.utils.build_watertap3 import watertap3_setup

warnings.filterwarnings('ignore')

__all__ = ['run_model',
            'run_and_return_model', 
            'run_model_no_print', 
            'watertap3_run', 
            'watertap3_binary_run',
            'case_study_constraints', 
            'get_ix_stash', 
            'fix_ix_stash',
            'print_ro_results', 
            'print_results', 
            'set_bounds',
            'set_flow', 'standard_ro_run',
            'get_ro_stash', 
            'fix_ro_stash', 
            'make_decision', 
            'connected_units', 
            'get_all_trains',  
            'run_all_trains',
            'run_watertap3_baseline', 
            'case_study_constraints'
            ]

def watertap3_run(m, solver='ipopt', return_df=False, tolerance=None, tee=False, \
                objective=True, flow_in=None, new_df_units=None, print_it=True, 
                only_system=False):
    
    '''
    Function to run WaterTAP3
    '''
    if print_it:
        print('\n=========================START WT3 MODEL RUN==========================')

    scenario = m.fs.train['scenario']
    case_study = m.fs.train['case_study']

    m = set_flow(m, flow_in)
    if print_it:
        run_model(m, solver=solver, objective=objective, tolerance=tolerance, tee=tee)
    else:
        run_model_no_print(m, solver=solver, objective=objective, tolerance=tolerance, tee=tee)
    m, ro_stash = standard_ro_run(m, objective=objective, solver=solver, print_it=print_it)

    if m.fs.has_ro:
        # m, ro_stash = get_ro_stash(m)

        m = watertap3_setup(case_study=case_study, scenario=scenario, 
                        ro_bounds=m.fs.ro_bounds, desired_recovery=m.fs.desired_recovery, 
                        new_df_units=new_df_units, print_it=print_it)
        m = get_case_study(m, new_df_units=new_df_units, print_it=print_it)
        m = set_flow(m, flow_in)
        if print_it:
            run_model(m, solver=solver, objective=objective, tolerance=tolerance, tee=tee)
        else:
            run_model_no_print(m, solver=solver, objective=objective, tolerance=tolerance, tee=tee)
            print(f'Degrees of Freedom = {degrees_of_freedom(m)}')
        m = fix_ro_stash(m, ro_stash)
        if m.fs.has_ix:
            m, ix_stash = get_ix_stash(m)
            m = fix_ix_stash(m, ix_stash)

    if print_it:
        run_model(m, solver=solver, objective=objective, \
            tolerance=tolerance, tee=tee, print_it=True, only_system=only_system)
    else:
        run_model_no_print(m, solver=solver, objective=objective, tolerance=tolerance, tee=tee)
        print(f'Degrees of Freedom = {degrees_of_freedom(m)}')

    m, df = get_results_table(m=m, case_study=case_study, scenario=scenario)
    if print_it:
        print('\n==========================END WT3 MODEL RUN===========================')

    if return_df:
        return m, df
    else:
        return m


def watertap3_binary_run(m, solver='gdpopt', return_df=False, objective=True,
                        flow_in=None, only_system=True, print_it=True):

    print('\n=========================START WT3 MODEL RUN==========================')

    scenario = m.fs.train['scenario']
    case_study = m.fs.train['case_study']

    all_trains_df = get_all_trains(m)
    m = set_flow(m, flow_in)
    try:
        run_model(m, solver=solver, objective=objective)
    except ValueError:
        print('BAD SOLVER STATUS -- Rebuilding and running with 1% higher flow rate...')
        flow_in_stash = flow_in
        flow_in = flow_in_stash * 1.01
        m = watertap3_setup(case_study=case_study, scenario=scenario, 
                        desired_recovery=m.fs.desired_recovery, ro_bounds=m.fs.ro_bounds, print_it=False)
        m = get_case_study(m)
        m = set_flow(m, flow_in)
        try:
            m = run_and_return_model(m, solver=solver, objective=objective, only_system=only_system, print_it=print_it)
        except ValueError:
            print('BAD SOLVER STATUS AGAIN -- Rebuilding and running with 5% higher flow rate...')
            flow_in = flow_in_stash * 1.05
            m = watertap3_setup(case_study=case_study, scenario=scenario, 
                        desired_recovery=m.fs.desired_recovery, ro_bounds=m.fs.ro_bounds, print_it=False)
            m = get_case_study(m, print_it=False)
            m = set_flow(m, flow_in)
            run_model(m, solver=solver, objective=objective)
    m, ro_stash = standard_ro_run(m, objective=objective, solver=solver, print_it=print_it)
    # m = fix_ro_stash(m, ro_stash)
    # run_model(m, solver=solver, objective=objective)
    # print_results(m, only_system=only_system)
    m = make_decision(m, case_study, scenario, flow_in=flow_in)
    # df_units = m.fs.df_units.copy()
    # m = watertap3_run(m, objective=objective, flow_in=flow_in, solver=solver, print_it=print_it, only_system=only_system, new_df_units=df_units)
    if m.fs.desired_recovery == 1:
        df_units = m.fs.df_units.copy()
        optimized_train = m.fs.optimized_train
        m = watertap3_run(m, objective=objective, flow_in=flow_in, solver=solver, \
            print_it=print_it, only_system=only_system, new_df_units=df_units)
        m.fs.all_trains_df = all_trains_df
        m.fs.optimized_train = optimized_train
    else:
        run_model(m, solver=solver, objective=objective)
        m, ro_stash = standard_ro_run(m, objective=objective, solver=solver)
        m = fix_ro_stash(m, ro_stash)
        run_model(m, solver=solver, objective=False)
    print_results(m, only_system=only_system)
    print(f'Optimized Train is Treatment Train {m.fs.optimized_train}\n')
    print(f'The following units were dropped:')
    optimized_train = m.fs.df_units.UnitName.to_list()
    all_trains = [t for t in m.fs.all_trains_df.binary_train.unique()]
    for t in all_trains:
        if t == m.fs.optimized_train:
            continue
        setattr(m.fs, f'df_units_train_{t}', m.fs.all_trains_df[m.fs.all_trains_df.binary_train == t].copy())
        train_df = getattr(m.fs, f'df_units_train_{t}')
        dropped_units = [u for u in train_df.UnitName.to_list() if u not in optimized_train]
        print(f'\n\tFrom Treatment Train {t}:')
        for u in dropped_units:
            print(f"\t\t{u.replace('_', ' ').swapcase()}")
    m = run_all_trains(m, flow_in=flow_in, only_system=only_system, objective=objective)
    m, df = get_results_table(m=m, case_study=case_study, scenario=scenario)

    print('\n==========================END WT3 MODEL RUN===========================')

    if return_df:
        return m, df
    else:
        return m        

def set_flow(m, flow_in):

    if isinstance(flow_in, dict):
        for source_name, flow in flow_in.items():
            source = getattr(m.fs, source_name)
            source.flow_vol_in.fix(flow)
    if isinstance(flow_in, float) or isinstance(flow_in, int):
        source_name = ast.literal_eval(m.fs.df_units.set_index(['Type']).loc['intake'].Parameter)['water_type'][0]
        source = getattr(m.fs, source_name)
        source.flow_vol_in.fix(flow_in)
    
    return m


def standard_ro_run(m, objective=True, solver='ipopt', print_it=True):

    m = set_bounds(m)
    if print_it:
        run_model(m, objective=objective, solver=solver)
    else:
        run_model_no_print(m, objective=objective, solver=solver)

    if m.fs.desired_recovery < 1:
        if m.fs.costing.system_recovery() > m.fs.desired_recovery:
            print('Running for desired recovery -->', m.fs.desired_recovery)
            m.fs.recovery_bound = Constraint(expr=
                m.fs.costing.system_recovery <= m.fs.desired_recovery)
            m.fs.recovery_bound1 = Constraint(expr=
                m.fs.costing.system_recovery >= m.fs.desired_recovery - 1.5)
            if print_it:
                run_model(m, objective=objective, solver=solver)
            else:
                run_model_no_print(m, objective=objective, solver=solver)
        else:
            print('System recovery already lower than desired recovery.'
                '\n\tDesired:', m.fs.desired_recovery, '\n\tCurrent:', m.fs.costing.system_recovery())
    
    m, ro_stash = get_ro_stash(m)
    return m, ro_stash

def run_all_trains(m, only_system=True, solver='ipopt', flow_in=None, objective=True,
                    incl_constituent_results=True):

    scenario = m.fs.train['scenario']
    case_study = m.fs.train['case_study']
    all_trains = [t for t in m.fs.all_trains_df.binary_train.unique()]
    m.fs.train_results = train_results = {}
    for train in all_trains:
        print(f'\nRunning train {train}...')
        train_results[train] = {}
        df_units = m.fs.all_trains_df[m.fs.all_trains_df.binary_train == train].copy()
        temp_m = watertap3_setup(case_study=case_study, scenario=scenario, 
                                new_df_units=df_units, print_it=False, ro_bounds=m.fs.ro_bounds,
                                desired_recovery=m.fs.desired_recovery)
        temp_m = get_case_study(temp_m, new_df_units=df_units, print_it=False)
        temp_m = watertap3_run(temp_m, objective=objective, flow_in=flow_in, \
                            new_df_units=df_units, solver=solver, print_it=False, only_system=True)
        # print(f'Degrees of Freedom = {degrees_of_freedom(temp_m)}')
        temp_scen = f'{case_study}_train_{train.swapcase()}'
        temp_m, df = get_results_table(m=temp_m, case_study=case_study, scenario=temp_scen,
                                        save=False, incl_constituent_results=incl_constituent_results)
        df['treatment_train'] = train
        train_results[train]['model'] = temp_m
        train_results[train]['results_df'] = df
        train_results[train]['LCOW'] = temp_m.fs.costing.LCOW()
        temp_m.fs.train['case_study'] = f'treatment_train_{train.swapcase()}'
        temp_m.fs.train['scenario'] = temp_scen
        print(f'***Train {train} results:***')
        print_results(temp_m, only_system=only_system)
    
    return m


def get_all_trains(m):

    trains = sorted([t for t in m.fs.df_units.binary_train.unique() if t != 'main'])
    all_trains_df = pd.DataFrame()
    for t in trains:

        train_df = m.fs.df_units[m.fs.df_units.binary_train.isin([t, 'main'])].copy()
        train_df['binary_train'] = t
        df = train_df.set_index(['UnitName']).copy()
        for unit in df.index.to_list():
            if 'split_fraction' not in str(df.loc[unit].Parameter):
                continue
            else:
                params = ast.literal_eval(df.loc[unit].Parameter)
                if all(x == 1 for x in params['split_fraction']):
                    series = df.loc[unit]
                    new_params = {k: v for k, v in params.items() if k != 'split_fraction'}
                    if new_params == {}: 
                        new_params = np.nan
                    series.Parameter = new_params
                    for to_unit, from_port in zip(series.ToUnitName.split(','), series.FromPort.split(',')):
                        if to_unit in df.index.to_list():
                            series.ToUnitName = to_unit
                            series.FromPort = from_port
                    df.loc[unit] = series
        all_trains_df = all_trains_df.append(df.reset_index())
    
    cols = ['binary_train', 'CaseStudy', 'Scenario', 'Unit', 'Type', 'UnitName',
                'ToUnitName', 'FromPort', 'Parameter', 'Reference']
    all_trains_df = all_trains_df[cols]
    return all_trains_df


def make_decision(m, case_study, scenario, flow_in=None):

    m.fs.units_to_drop = units_to_drop = []
    m.fs.units_to_keep = units_to_keep = []
    all_dropped_units = []

    all_trains_df = get_all_trains(m)

    for splitter_name, splitter_info in m.fs.all_splitters.items():
        remove_units = []
        if splitter_info['indicator']:
            splitter = getattr(m.fs, splitter_name)
            from_unit = splitter._split_from_unit
            for out in splitter.outlet_list:
                disjunct = getattr(splitter, f'disjunct_{out}')
                outlet = getattr(splitter, out)
                if bool(disjunct.indicator_var):
                    chosen_unit = outlet.to_unit
                else:
                    remove_units.append(outlet.to_unit)
                    all_dropped_units.append(outlet.to_unit)
        else:
            continue

        df_units = m.fs.df_units.set_index(['UnitName']).drop(index=remove_units).copy()
        from_unit_series = df_units.loc[from_unit].copy()
        from_unit_params = ast.literal_eval(from_unit_series.Parameter)
        from_unit_series.Parameter = str({k: v for k, v in from_unit_params.items() if k != 'split_fraction'})
        to_unit_name = []
        from_port = []
        for unit, port in zip(from_unit_series.ToUnitName.split(','), from_unit_series.FromPort.split(',')):
            if unit == chosen_unit:
                to_unit_name.append(unit)
                from_port.append(port)
                continue
            if port == 'waste':
                to_unit_name.append(unit)
                from_port.append(port)
                continue
        from_unit_series.ToUnitName = ','.join(to_unit_name)
        from_unit_series.FromPort = ','.join(from_port)
        df_units.loc[from_unit] = from_unit_series
        df_units.reset_index(inplace=True)

        temp_pfd_dict = get_pfd_dict(df_units)

        for remove in remove_units: 
            start_u = m.fs.pfd_dict[remove]['ToUnitName']
            if isinstance(start_u, list):
                temp_drop = connected_units(start_u, temp_pfd_dict, units=[s for s in start_u])
            else:
                temp_drop = connected_units(start_u, temp_pfd_dict, units=[start_u])
            units_to_drop += temp_drop
        temp_keep = connected_units(from_unit, temp_pfd_dict, units=[])
        units_to_keep += temp_keep
    
    units_to_keep = list(set(units_to_keep))
    units_to_drop = [u for u in list(set(units_to_drop)) if u not in units_to_keep]
    all_dropped_units += units_to_drop
    df_units = df_units.set_index('UnitName').drop(index=units_to_drop).copy()
    df_units.reset_index(inplace=True)

    optimized_train = df_units[df_units.binary_train != 'main'].binary_train.unique()[0]
    
    m = watertap3_setup(case_study=case_study, scenario=scenario, 
                    new_df_units=df_units, print_it=False, ro_bounds=m.fs.ro_bounds,
                    desired_recovery=m.fs.desired_recovery)
    m.fs.optimized_train = optimized_train
    m.fs.all_dropped_units = all_dropped_units
    m.fs.all_trains_df = all_trains_df
    m = get_case_study(m, new_df_units=df_units, print_it=False)
    m = set_flow(m, flow_in)


    return m

def run_model(m, solver='ipopt', tolerance=None, tee=False, objective=False, 
                max_attempts=3, print_it=False, initial_run=True, mip_solver='glpk', only_system=False):
    
    '''
    Function used to attempt model solve.
    '''

    if initial_run:
        financials.get_system_costing(m.fs)

    TransformationFactory('network.expand_arcs').apply_to(m)

    if objective:
        if objective in ['TCI', 'elec', 'fixed_op', 'chem', 'other_onm', 'total_op']:
            obj = getattr(m.fs.costing, f'LCOW_{objective}')
        elif objective == 'electricity_intensity':
            obj = (m.fs.costing.electricity_cost_annual / m.fs.costing_param.electricity_price) / \
                    m.fs.costing.treated_water
        else:
            obj = m.fs.costing.LCOW
        m.fs.objective_function = Objective(expr=obj)

    model_solver = SolverFactory(solver)
    if tolerance and solver == 'ipopt':
        model_solver.options['tol'] = tolerance

    logging.getLogger('pyomo.core').setLevel(logging.ERROR)

    print('.................................')
    print('\nDegrees of Freedom:', degrees_of_freedom(m))
    if solver == 'gdpopt':
        m.fs.results = results = model_solver.solve(m, tee=tee, mip_solver=mip_solver)
    else:
        m.fs.results = results = model_solver.solve(m, tee=tee)
    print(f'\nInitial solve attempt {results.solver.termination_condition.swapcase()}')

    attempt_number = 1
    while ((m.fs.results.solver.termination_condition in \
        ['infeasible', 'maxIterations', 'unbounded', 'other']) & (attempt_number <= max_attempts)):
        print(f'\nAttempt {attempt_number}:')
        if solver == 'gdpopt':
            m.fs.results = results = model_solver.solve(m, tee=tee, mip_solver=mip_solver)
        else:
            m.fs.results = results = model_solver.solve(m, tee=tee)
        print(f'\n\tWaterTAP3 solver returned {results.solver.termination_condition.swapcase()} solution...')
        # if results.solver.termination_condition == 'infeasible':
        #     print(infeas.log_infeasible_bounds(m))
        #     print(infeas.log_infeasible_constraints(m))
        attempt_number += 1
    if m.fs.results.solver.termination_condition != 'optimal':
        raise Exception(f'\n\tMODEL RUN ABORTED:'
              f'\n\tWT3 solution is {m.fs.results.solver.termination_condition.swapcase()}'
              f'\n\tModel did not solve optimally after 3 attempts. No results are saved.'
              f'\n\tCheck model setup and initial conditions and retry.')
    print(f'\nWaterTAP3 solution {results.solver.termination_condition.swapcase()}\n')
    print('.................................')

    if print_it:
        print_results(m, only_system=only_system)


def run_and_return_model(m, solver='ipopt', tolerance=None, tee=False, objective=False, 
                        max_attempts=3, print_it=False, initial_run=True, mip_solver='glpk', only_system=False):

    '''
    Function to attempt model solve and return model object.
    '''

    if initial_run:
        financials.get_system_costing(m.fs)

    TransformationFactory('network.expand_arcs').apply_to(m)

    if objective:
        if objective in ['TCI', 'elec', 'fixed_op', 'chem', 'other_onm', 'total_op']:
            obj = getattr(m.fs.costing, f'LCOW_{objective}')
        elif objective == 'electricity_intensity':
            obj = (m.fs.costing.electricity_cost_annual / m.fs.costing_param.electricity_price) / \
                    m.fs.costing.treated_water
        else:
            obj = m.fs.costing.LCOW
        m.fs.objective_function = Objective(expr=obj)

    model_solver = SolverFactory(solver)
    if tolerance and solver == 'ipopt':
        model_solver.options['tol'] = tolerance

    logging.getLogger('pyomo.core').setLevel(logging.ERROR)

    print('.................................')
    print('\nDegrees of Freedom:', degrees_of_freedom(m))
    if solver == 'gdpopt':
        m.fs.results = results = model_solver.solve(m, tee=tee, mip_solver=mip_solver)
    else:
        m.fs.results = results = model_solver.solve(m, tee=tee)
    
    print(f'\nInitial solve attempt {results.solver.termination_condition.swapcase()}')

    attempt_number = 1
    while ((m.fs.results.solver.termination_condition in \
        ['infeasible', 'maxIterations', 'unbounded', 'other']) & (attempt_number <= max_attempts)):
        print(f'\nAttempt {attempt_number}:')
        if solver == 'gdpopt':
            m.fs.results = results = model_solver.solve(m, tee=tee, mip_solver=mip_solver)
        else:
            m.fs.results = results = model_solver.solve(m, tee=tee)
        print(f'\n\tWaterTAP3 solver returned {results.solver.termination_condition.swapcase()} solution...')
        attempt_number += 1
    if m.fs.results.solver.termination_condition != 'optimal':
        raise Exception(f'\n\tMODEL RUN ABORTED:'
              f'\n\tWT3 solution is {m.fs.results.solver.termination_condition.swapcase()}'
              f'\n\tModel did not solve optimally after 3 attempts. No results are saved.'
              f'\n\tCheck model setup and initial conditions and retry.')
    print(f'\nWaterTAP3 solution {results.solver.termination_condition.swapcase()}\n')
    print('.................................')

    if print_it:
        print_results(m, only_system=only_system)

    return m


def run_model_no_print(m, solver='ipopt', tolerance=None, tee=False, objective=False, 
                        max_attempts=3, initial_run=True, mip_solver='glpk', return_model=False):

    if initial_run:
        financials.get_system_costing(m.fs)

    TransformationFactory('network.expand_arcs').apply_to(m)

    if objective:
        if objective in ['TCI', 'elec', 'fixed_op', 'chem', 'other_onm', 'total_op']:
            obj = getattr(m.fs.costing, f'LCOW_{objective}')
        elif objective == 'electricity_intensity':
            obj = (m.fs.costing.electricity_cost_annual / m.fs.costing_param.electricity_price) / \
                    m.fs.costing.treated_water
        else:
            obj = m.fs.costing.LCOW
        m.fs.objective_function = Objective(expr=obj)

    model_solver = SolverFactory(solver)
    if tolerance and solver == 'ipopt':
        model_solver.options['tol'] = tolerance
    # m.fs.solver = solver = SolverFactory('glpk')

    logging.getLogger('pyomo.core').setLevel(logging.ERROR)

    if solver == 'gdpopt':
        m.fs.results = model_solver.solve(m, tee=tee, mip_solver=mip_solver)
    else:
        m.fs.results = model_solver.solve(m, tee=tee)

    attempt_number = 1
    while ((m.fs.results.solver.termination_condition in \
        ['infeasible', 'maxIterations', 'unbounded']) & (attempt_number <= max_attempts)):
        if solver == 'gdpopt':
            m.fs.results = model_solver.solve(m, tee=tee, mip_solver=mip_solver)
        else:
            m.fs.results = model_solver.solve(m, tee=tee)
        attempt_number += 1
    if m.fs.results.solver.termination_condition != 'optimal':
        raise Exception(f'\n\tMODEL RUN ABORTED:'
              f'\n\tWT3 solution is {m.fs.results.solver.termination_condition.swapcase()}'
              f'\n\tModel did not solve optimally after 3 attempts. No results are saved.'
              f'\n\tCheck model setup and initial conditions and retry.')
    if return_model:
        return m


def get_ix_stash(m):
    m.fs.ix_stash = ix_stash = {}
    df = m.fs.df_units.set_index(['UnitName'])
    for u in df.index:
        unit_module = df.loc[u].Unit
        if unit_module == 'ion_exchange':
            unit = getattr(m.fs, u)
            ix_stash[u] = {
                        'sfr': unit.sfr(),
                        'resin_depth': unit.resin_depth(),
                        'column_diam': unit.column_diam(),
                        'num_columns': unit.num_columns()
                        }

    return m, ix_stash


def fix_ix_stash(m, ix_stash, only_num_cols=False):
    if only_num_cols:
        for ix in ix_stash.keys():
            unit = getattr(m.fs, ix)
            unit.num_columns.fix(ix_stash[ix]['num_columns'])
        return m
    else:
        for ix in ix_stash.keys():
            unit = getattr(m.fs, ix)
            unit.sfr.fix(ix_stash[ix]['sfr'])
            # unit.num_columns.fix(ix_stash[ix]['num_columns'])
            unit.resin_depth.fix(ix_stash[ix]['resin_depth'])
            # unit.column_diam.fix(ix_stash[ix]['column_diam'])

        return m


def get_ro_stash(m):
    m.fs.ro_stash = ro_stash = {}
    for k, v in m.fs.pfd_dict.items():
        if v['Unit'] == 'reverse_osmosis':
            unit = getattr(m.fs, k)
            ro_stash[k] = {
                    'feed.pressure': unit.feed.pressure[0](),
                    'membrane_area': unit.membrane_area[0](),
                    'a': unit.a[0](),
                    'b': unit.b[0]()
                    }
    return m, ro_stash


def fix_ro_stash(m, ro_stash):
    for ro in ro_stash.keys():
        unit = getattr(m.fs, ro)
        unit.feed.pressure.fix(ro_stash[ro]['feed.pressure'])
        unit.membrane_area.fix(ro_stash[ro]['membrane_area'])
        unit.a.fix(ro_stash[ro]['a'])
        unit.b.fix(ro_stash[ro]['b'])

    return m


def print_results(m, only_system=False):
    case_study = m.fs.train['case_study']
    scenario = m.fs.train['scenario']
    case_study_print = case_study.replace('_', ' ').swapcase()
    scenario_print = scenario.replace('_', ' ').swapcase()
    print(f'\n{case_study_print}: {scenario_print}')
    print('=========================SYSTEM LEVEL RESULTS=========================')
    print('LCOW ($/m3):', round(value(m.fs.costing.LCOW()), 5))
    print('Total Capital Investment ($MM):', round(value(m.fs.costing.capital_investment_total()), 3))
    print('Total Annual Operating Costs ($MM/yr):', round(value(m.fs.costing.operating_cost_annual()), 3))
    print('Annual Fixed Operating Cost ($MM/yr):', round(value(m.fs.costing.fixed_op_cost_annual()), 3))
    print('Annual Chemicals Cost ($MM/yr):', round(value(m.fs.costing.cat_and_chem_cost_annual()), 3))
    print('Annual Electricity Costs ($MM/yr):', round(value(m.fs.costing.electricity_cost_annual()), 3))
    print('Annual Other Variable Costs ($MM/yr):', round(value(m.fs.costing.other_var_cost_annual()), 3))
    print('Total flow in (m3/s):', round(value(m.fs.sys_flow_in), 3))
    print('Total flow out (m3/s):', round(value(m.fs.costing.treated_water()), 3))
    print('Total flow waste (m3/s):', round(value(m.fs.sys_flow_in), 3) - round(value(m.fs.costing.treated_water()), 3))
    print('Total water recovery (%):', round(value(100 * m.fs.costing.system_recovery()), 3))
    print('Electricity intensity (kWh/m3):', round(value(m.fs.costing.electricity_intensity()), 3))
    print('Electricity portion of LCOW (%):', round(value(100 * m.fs.costing.elec_frac_LCOW()), 3))
    print('======================================================================')
    if only_system:
        return
    print('\n=========================UNIT PROCESS RESULTS=========================\n')
    for unit in m.fs.df_units.UnitName:
        b_unit = getattr(m.fs, unit)
        print(f'\n{b_unit.unit_pretty_name}:')
        print('\tTotal Capital Investment ($MM):', round(value(b_unit.costing.total_cap_investment()), 5))
        print('\tAnnual O&M ($MM/yr):', round(value(b_unit.costing.annual_op_main_cost), 5))
        print('\tAnnual Fixed O&M ($MM/yr):', round(value(b_unit.costing.total_fixed_op_cost), 5))
        print('\tAnnual Chemical Cost ($MM/yr):', round(value(b_unit.costing.cat_and_chem_cost), 5))
        print('\tAnnual Other O&M Cost ($MM/yr):', round(value(b_unit.costing.other_var_cost), 5))
        print('\tAnnual Electricity Cost ($MM/yr):', round(value(b_unit.costing.electricity_cost), 5))
        print('\tElectricity Intensity (kWh/m3):', round(value(b_unit.costing.electricity_intensity()), 5))
        print('\tUnit LCOW ($/m3):', round(value(b_unit.LCOW()), 5))
        print('\tFlow In (m3/s):', round(value(b_unit.flow_vol_in[0]()), 5))
        print('\tFlow Out (m3/s):', round(value(b_unit.flow_vol_out[0]()), 5))
        if hasattr(b_unit, 'flow_vol_waste'):
            print('\tFlow Waste (m3/s):', round(value(b_unit.flow_vol_waste[0]()), 5))
        if b_unit.unit_type == 'reverse_osmosis':
            print_ro_results(m, b_unit.unit_name)
        if b_unit.unit_type == 'chlorination':
            print('\tChlorine dose (mg/L):', round(b_unit.dose, 3))
        if b_unit.unit_type in ['brine_concentrator', 'evaporation_pond', 'landfill', 'landfill_zld', \
            'microfiltration', 'ultrafiltration', 'nanofiltration'] or b_unit.unit_kind == 'use':
            try:
                print('\tTDS in (mg/L):', round(value(b_unit.conc_mass_in[0, 'tds']) * 1000, 1))
                print('\tTDS out (mg/L):', round(value(b_unit.conc_mass_out[0, 'tds']) * 1000, 1))
                if hasattr(b_unit, 'conc_mass_waste'):
                    print('\tTDS waste (mg/L):', round(value(b_unit.conc_mass_waste[0, 'tds']) * 1000, 1))
            except:
                print(f'\tNO TDS INTO {b_unit.unit_pretty_name}')
            if b_unit.unit_type == 'evaporation_pond':
                print(f'\tPond Area (acres): {round(b_unit.area[0](), 3)}')
            print('\tWater Recovery (%):', round(value((b_unit.flow_vol_out[0]() / b_unit.flow_vol_in[0]())), 5) * 100)
            continue
        if b_unit.unit_type == 'uv_aop':
            print(f'\tUV Dose (mJ/cm2): {round(b_unit.uv_dose, 1)}')
            print(f'\tUVT (%): {round(b_unit.uvt_in, 3) * 100}%')
            if hasattr(b_unit, 'water_recovery'):
                print('\tWater Recovery (%):', round(value(b_unit.water_recovery[0]()), 5) * 100)
            continue
        elif b_unit.unit_type != 'reverse_osmosis':
            if hasattr(b_unit, 'water_recovery'):
                print('\tWater Recovery (%):', round(value(b_unit.water_recovery[0]()), 5) * 100)

    print('\n======================================================================\n')


def print_ro_results(m, ro_name):
    pressures = []
    recovs = []
    areas = []
    num_mems = []
    kws = []
    kss = []
    fluxs = []
    flow_ins = []
    flow_outs = []
    tds_in = []
    tds_out = []
    tds_waste = []

    unit = getattr(m.fs, ro_name)
    pressures.append(unit.feed.pressure[0]())
    recovs.append(unit.ro_recovery())
    areas.append(unit.membrane_area[0]())
    num_mems.append(unit.num_membranes())
    kws.append(unit.a[0]())
    kss.append(unit.b[0]())
    fluxs.append(unit.flux_lmh)
    flow_ins.append(unit.flow_vol_in[0]())
    flow_outs.append(unit.flow_vol_out[0]())
    tds_in.append(unit.conc_mass_in[0, 'tds']())
    tds_out.append(unit.conc_mass_out[0, 'tds']())
    tds_waste.append(unit.conc_mass_waste[0, 'tds']())
    # print(f'.. {print_name}:')
    print(f'\n\t.....{unit.unit_pretty_name} OPERATIONAL PARAMETERS.....')
    print(f'\tPressure = {round(pressures[-1], 2)} bar = {round(pressures[-1] * 14.5038)} psi')
    print(f'\tArea = {round(areas[-1])} m2 ---> {round(num_mems[-1])} membrane modules')
    print(f'\tFlux = {round(value(fluxs[-1]), 1)} LMH')
    print(f'\tTDS in = {round(value(tds_in[-1]) * 1000, 3)} mg/L')
    print(f'\tTDS out = {round(value(tds_out[-1]) * 1000, 3)} mg/L')
    print(f'\tTDS waste = {round(value(tds_waste[-1]) * 1000, 3)} mg/L')
    print(f'\tFlow in = {round(unit.flow_vol_in[0](), 5)} m3/s = '
          f'{round(pyunits.convert(unit.flow_vol_in[0], to_units=pyunits.Mgallons / pyunits.day)(), 5)} MGD = '
          f'{round(pyunits.convert(unit.flow_vol_in[0], to_units=pyunits.gallons / pyunits.min)(), 2)} gpm')
    print(f'\tFlow out = {round(unit.flow_vol_out[0](), 5)} m3/s = '
          f'{round(pyunits.convert(unit.flow_vol_out[0], to_units=pyunits.Mgallons / pyunits.day)(), 5)} MGD = '
          f'{round(pyunits.convert(unit.flow_vol_out[0], to_units=pyunits.gallons / pyunits.min)(), 2)} gpm')
    print(f'\tFlow waste = {round(unit.flow_vol_waste[0](), 5)} m3/s = '
          f'{round(pyunits.convert(unit.flow_vol_waste[0], to_units=pyunits.Mgallon / pyunits.day)(), 5)} MGD = '
          f'{round(pyunits.convert(unit.flow_vol_waste[0], to_units=pyunits.gallons / pyunits.min)(), 2)} gpm')
    print(f'\tWater Permeability = {kws[-1]} m/(bar.hr)')
    print(f'\tSalt Permeability = {kss[-1]} m/hr')
    print(f'\tRO Recovery = {round(recovs[-1], 3) * 100}%')

def set_bounds(m):

    if m.fs.ro_bounds == 'swro':
        feed_flux_max = 45  # lmh
        feed_flux_min = 10  # lmh
        a = [2, 7]
        b = [0.2, 0.7]
        max_pressure = 85
        min_area = 500
        min_pressure = 5
    # if m.fs.ro_bounds == 'hpro':
    #     feed_flux_max = 100  # lmh
    #     feed_flux_min = 1  # lmh
    #     a = [2, 9]
    #     b = [0.2, 0.9]
    #     max_pressure = 120
    #     min_area = 500
    #     min_pressure = 5
    # if m.fs.ro_bounds == 'bwro':
    #     feed_flux_max = 30  # lmh
    #     feed_flux_min = 8  # lmh
    #     a = [2, 7]
    #     b = [0.2, 0.7]
    #     max_pressure = 15
    #     min_area = 500
    #     min_pressure = 5
    else:
        feed_flux_max = 30  # lmh
        feed_flux_min = 8  # lmh
        a = [2, 7]
        b = [0.2, 0.7]
        max_pressure = 25
        min_area = 250
        min_pressure = 5

    q = 1
    for key in m.fs.pfd_dict.keys():
        if m.fs.pfd_dict[key]['Unit'] == 'reverse_osmosis':

            setattr(m, ('flux_constraint%s' % q), Constraint(
                    expr=getattr(m.fs, key).pure_water_flux[0] * 3600 <= feed_flux_max))
            q += 1
            setattr(m, ('flux_constraint%s' % q), Constraint(
                    expr=getattr(m.fs, key).pure_water_flux[0] * 3600 >= feed_flux_min))
            q += 1
            setattr(m, ('flux_constraint%s' % q), Constraint(
                    expr=getattr(m.fs, key).membrane_area[0] >= min_area))
            q += 1
            setattr(m, ('flux_constraint%s' % q), Constraint(
                    expr=getattr(m.fs, key).feed.pressure[0] <= max_pressure))

            # q += 1
            # setattr(m, ('flux_constraint%s' % q), Constraint(
            #         expr=getattr(m.fs, key).feed.pressure[0] >= min_pressure))
            q += 1
            setattr(m.fs, ('flux_constraint%s' % q), Constraint(
                    expr=getattr(m.fs, key).a[0] <= a[1]))
            q += 1
            setattr(m.fs, ('flux_constraint%s' % q), Constraint(
                    expr=getattr(m.fs, key).a[0] >= a[0]))
            q += 1
            setattr(m.fs, ('flux_constraint%s' % q), Constraint(
                    expr=getattr(m.fs, key).b[0] <= b[1]))
            q += 1
            setattr(m.fs, ('flux_constraint%s' % q), Constraint(
                    expr=getattr(m.fs, key).b[0] >= b[0]))
            q += 1

    return m



def connected_units(u, temp_pfd_dict, units=[]):
    if isinstance(u, list):
        for unit in u:
            if temp_pfd_dict[unit]['binary_train'] == 'main':
                continue
            next_unit = temp_pfd_dict[unit]['ToUnitName']
            if next_unit in units or next_unit is np.nan:
                continue
            if isinstance(next_unit, list):
                units += [n for n in next_unit]
            else:
                units.append(next_unit)
            connected_units(next_unit, temp_pfd_dict, units=units)
        return units
    if u is np.nan:
        return units
    # try:
    next_unit = temp_pfd_dict[u]['ToUnitName']
    # except KeyError as e:
    #     units.append(u)
    #     return units
    if next_unit is np.nan:
        # print('next_unit is np.nan')
        return units
    if isinstance(next_unit, list):
        for unit in next_unit:
            if unit in units:
                continue
            units.append(unit)
            connected_units(unit, temp_pfd_dict, units=units)
        return units
    else:
        if next_unit in units:
            return units
        units.append(next_unit)
        connected_units(next_unit, temp_pfd_dict, units=units)
        return units



def run_watertap3_baseline(m, solver='ipopt', 
                    return_df=False, tolerance=None, tee=False, objective=True):
    
    '''
    Function to run old WaterTAP3 Baseline case studies
    '''
    print('\n=========================START WT3 MODEL RUN==========================')
    scenario = m.fs.train['scenario']
    case_study = m.fs.train['case_study']
    reference = m.fs.train['reference']

    run_model(m, solver=solver, objective=objective, tolerance=tolerance, tee=tee)

    if m.fs.choose:
        print(f'********TREATMENT TRAIN CONTAINS DECISION VARIABLE********')
        print('Removing non-optimal unit processes...\n\n')
        all_trains_df = get_all_trains(m)
        m = make_decision(m, case_study, scenario)
        print('The following units were dropped:')
        for dropped_unit in m.fs.all_dropped_units:
            print(f"\t{dropped_unit.replace('_', ' ').swapcase()}")
        print('\n=======================OPTIMIZED TREATMENT TRAIN=======================')
        run_model(m, solver=solver, objective=objective, tolerance=tolerance, tee=tee)

    if m.fs.has_ix:
        m, ix_stash = get_ix_stash(m)
        print('Initial IX solve OK...\nFixing number IX columns...')
        m = fix_ix_stash(m, ix_stash)
        run_model(m, solver=solver, objective=objective, tolerance=tolerance)

    m = case_study_constraints(m, case_study, scenario)

    if m.fs.has_ro:
        if case_study == 'upw':
            m.fs.splitter2.split_fraction_constr = Constraint(expr=
                sum(m.fs.splitter2.split_fraction_vars) <= 1.001)
            m.fs.splitter2.split_fraction_constr2 = Constraint(expr=
                sum(m.fs.splitter2.split_fraction_vars) >= 0.999)
        m = set_bounds(m)
        run_model(m, solver=solver, objective=objective, tolerance=tolerance)

    if m.fs.desired_recovery < 1:
        if m.fs.costing.system_recovery() > m.fs.desired_recovery:
            print('Running for desired recovery -->', m.fs.desired_recovery)
            m.fs.recov_ub = Constraint(expr=
                m.fs.costing.system_recovery <= m.fs.desired_recovery)
            m.fs.recov_lb = Constraint(expr=
                m.fs.costing.system_recovery >= m.fs.desired_recovery - 1.5)
            run_model(m, objective=objective, tolerance=tolerance)
        
        else:
            print('System recovery already lower than desired recovery.'
                  '\n\tDesired:', m.fs.desired_recovery, '\n\tCurrent:', m.fs.costing.system_recovery())

    if case_study == 'uranium':
        ur_list = []
        ur_list.append(m.fs.ion_exchange.removal_fraction[0, 'tds']())
        ur_list.append(m.fs.ion_exchange.anion_res_capacity[0]())
        ur_list.append(m.fs.ion_exchange.cation_res_capacity[0]())

    if case_study == 'upw':
        m.fs.upw_list = upw_list = []
        upw_list.append(m.fs.splitter2.split_fraction_outlet_1[0]())
        upw_list.append(m.fs.splitter2.split_fraction_outlet_2[0]())

    if m.fs.has_ro:
        m, ro_stash = get_ro_stash(m)
        ###### RESET BOUNDS AND DOUBLE CHECK RUN IS OK SO CAN GO INTO SENSITIVITY #####
        if m.fs.new_case_study:
            new_df_units = m.fs.df_units.copy()
            all_dropped_units = m.fs.all_dropped_units
            m = watertap3_setup(dynamic=False, case_study=case_study, scenario=scenario, 
                                new_df_units=new_df_units, ro_bounds=m.fs.ro_bounds, 
                                desired_recovery=m.fs.desired_recovery)
            m.fs.all_dropped_units = all_dropped_units
            m = get_case_study(m, new_df_units=new_df_units)


        else:
            m = watertap3_setup(dynamic=False, case_study=case_study, scenario=scenario, 
                                ro_bounds=m.fs.ro_bounds, desired_recovery=m.fs.desired_recovery)
            m = get_case_study(m)


        if case_study == 'gila_river' and scenario != 'baseline':
            m.fs.evaporation_pond.water_recovery.fix(0.895)

        if case_study == 'upw':
            run_model(m, solver=solver, objective=objective, tolerance=tolerance)
            m.fs.upw_list = upw_list
            m.fs.media_filtration.water_recovery.fix(0.9)
            m.fs.splitter2.split_fraction_outlet_1.fix(upw_list[0])
            m.fs.splitter2.split_fraction_outlet_2.fix(upw_list[1])

        if case_study == 'ocwd':  
            m.fs.microfiltration.water_recovery.fix(0.9)

        if case_study == 'uranium':
            m.fs.ion_exchange.removal_fraction[0, 'tds'].fix(ur_list[0])
            m.fs.ion_exchange.anion_res_capacity.fix(ur_list[1])
            m.fs.ion_exchange.cation_res_capacity.fix(ur_list[2])

        if case_study == 'irwin':
            run_model(m, solver=solver, objective=objective, tolerance=tolerance)
            m.fs.brine_concentrator.water_recovery.fix(0.8)
        
        run_model(m, solver=solver, objective=objective, tolerance=tolerance)
        m = fix_ro_stash(m, ro_stash)
        m.fs.objective_function.deactivate()
        # m = fix_ro_stash(m, ro_stash)
        if m.fs.has_ix:
        #     m, ix_stash = get_ix_stash(m)
            m = fix_ix_stash(m, ix_stash)

    run_model(m, solver=solver, objective=False, print_it=True, tolerance=tolerance)

    m, df = get_results_table(m=m, case_study=case_study, scenario=scenario)

    print('\n==========================END WT3 MODEL RUN===========================')

    if return_df:
        return m, df
    else:
        return m


def case_study_constraints(m, case_study, scenario):
    if case_study == 'upw':
        m.fs.media_filtration.water_recovery.fix(0.9)
        m.fs.reverse_osmosis.eq1_upw = Constraint(expr=
            m.fs.reverse_osmosis.flow_vol_out[0] <= 0.05678 * 1.01)
        m.fs.reverse_osmosis.eq2_upw = Constraint(expr=
            m.fs.reverse_osmosis.flow_vol_out[0] >= 0.05678 * 0.99)
        m.fs.reverse_osmosis.eq3_upw = Constraint(expr=
            m.fs.reverse_osmosis.flow_vol_waste[0] <= 0.04416 * 1.01)
        m.fs.reverse_osmosis.eq4_upw = Constraint(expr=
            m.fs.reverse_osmosis.flow_vol_waste[0] >= 0.04416 * 0.99)
        m.fs.reverse_osmosis_2.eq1 = Constraint(expr=
            m.fs.reverse_osmosis_2.flow_vol_out[0] <= 0.01262 * 1.01)
        m.fs.reverse_osmosis_2.eq2 = Constraint(expr=
            m.fs.reverse_osmosis_2.flow_vol_out[0] >= 0.01262 * 0.99)
        m.fs.ro_stage.eq1_upw = Constraint(expr=
            m.fs.ro_stage.flow_vol_out[0] <= 0.03154 * 1.01)
        m.fs.ro_stage.eq2_upw = Constraint(expr=
            m.fs.ro_stage.flow_vol_out[0] >= 0.03154 * 0.99)

        if scenario not in ['baseline']:
            m.fs.to_zld_constr1 = Constraint(expr=
                m.fs.to_zld.flow_vol_in[0] >= 0.99 * 0.0378)
            m.fs.to_zld_constr2 = Constraint(expr=
                m.fs.to_zld.flow_vol_in[0] <= 1.01 * 0.0378)

    if case_study == 'uranium':
        m.fs.ro_production.eq1 = Constraint(expr=
            m.fs.ro_production.flow_vol_out[0] <= 
            (0.7 * m.fs.ro_production.flow_vol_in[0]) * 1.01)
        m.fs.ro_production.eq2 = Constraint(expr=
            m.fs.ro_production.flow_vol_out[0] >= 
            (0.7 * m.fs.ro_production.flow_vol_in[0]) * 0.99)
        m.fs.ro_restore_stage.eq3 = Constraint(expr=
            m.fs.ro_restore_stage.flow_vol_out[0] <= 
            (0.5 * m.fs.ro_restore_stage.flow_vol_in[0]) * 1.01)
        m.fs.ro_restore_stage.eq4 = Constraint(expr=
            m.fs.ro_restore_stage.flow_vol_out[0] >= 
            (0.5 * m.fs.ro_restore_stage.flow_vol_in[0]) * 0.99)
        m.fs.ro_restore.eq5 = Constraint(expr=
            m.fs.ro_restore.flow_vol_out[0] <= 
            (0.75 * m.fs.ro_restore.flow_vol_in[0]) * 1.01)
        m.fs.ro_restore.eq6 = Constraint(expr=
            m.fs.ro_restore.flow_vol_out[0] >= 
            (0.75 * m.fs.ro_restore.flow_vol_in[0]) * 0.99)

    if case_study == 'gila_river':
        if 'reverse_osmosis' in m.fs.pfd_dict.keys():
            m.fs.reverse_osmosis.recov1 = Constraint(expr=
                m.fs.reverse_osmosis.flow_vol_out[0] <= 
                (0.59 * m.fs.reverse_osmosis.flow_vol_in[0]) * 1.01)
            m.fs.reverse_osmosis.recov2 = Constraint(expr=
                m.fs.reverse_osmosis.flow_vol_out[0] >= 
                (0.59 * m.fs.reverse_osmosis.flow_vol_in[0]) * 0.99)

        if scenario != 'baseline':
            m.fs.evaporation_pond.water_recovery.fix(0.895)

    if case_study == 'cherokee':

        if 'boiler_ro' in m.fs.pfd_dict.keys():
            m.fs.boiler_ro.recov1 = Constraint(expr=
                m.fs.boiler_ro.flow_vol_out[0] <= \
                    (0.75 * m.fs.boiler_ro.flow_vol_in[0]) * 1.01)
            m.fs.boiler_ro.recov2 = Constraint(expr=
                m.fs.boiler_ro.flow_vol_out[0] >= \
                    (0.75 * m.fs.boiler_ro.flow_vol_in[0]) * 0.99)

        if 'reverse_osmosis_a' in m.fs.pfd_dict.keys():
            m.fs.reverse_osmosis_a.recov1 = Constraint(expr=
                m.fs.reverse_osmosis_a.flow_vol_out[0] >= \
                    (0.95 * m.fs.reverse_osmosis_a.flow_vol_in[0]))

    if case_study == 'san_luis':
        if scenario in ['baseline', 'dwi', '1p5_mgd', '3_mgd', '5_mgd']:
            m.fs.reverse_osmosis_1.feed.pressure.fix(25.5)
            m.fs.reverse_osmosis_2.feed.pressure.fix(36)

    if case_study == 'kbhdp':
        m.fs.ro_recovery_constr1 = Constraint(expr=
            (m.fs.ro_first_stage.flow_vol_out[0] + 
                m.fs.ro_second_stage.flow_vol_out[0]) /
                m.fs.ro_first_stage.flow_vol_in[0] <= 0.83)
        m.fs.ro_recovery_constr2 = Constraint(expr=
            (m.fs.ro_first_stage.flow_vol_out[0] + 
                m.fs.ro_second_stage.flow_vol_out[0]) / 
                m.fs.ro_first_stage.flow_vol_in[0] >= 0.81)
        m.fs.ro1_press_constr2 = Constraint(expr=
            m.fs.ro_first_stage.feed.pressure[0] <= 14)
        m.fs.ro2_press_constr2 = Constraint(expr=
            m.fs.ro_second_stage.feed.pressure[0] <= 18)
        m.fs.ro_area_constr1 = Constraint(expr=
            m.fs.ro_first_stage.membrane_area[0] /
            m.fs.ro_second_stage.membrane_area[0] <= 2.1)
        m.fs.ro_area_constr2 = Constraint(expr=
            m.fs.ro_first_stage.membrane_area[0] / 
            m.fs.ro_second_stage.membrane_area[0] >= 1.9)  

    if case_study == 'emwd':
        if scenario in ['baseline', 'dwi']:
            m.fs.manifee_area_constr = Constraint(expr=
                m.fs.menifee_a.membrane_area[0] == 
                m.fs.menifee_b.membrane_area[0])
            m.fs.perris_area_constr = Constraint(expr=
                m.fs.perris_i_a.membrane_area[0] == 
                m.fs.perris_i_b.membrane_area[0])
            m.fs.menifee_pressure_constr1 = Constraint(expr=
                m.fs.menifee_a.feed.pressure[0] <= 14)
            m.fs.menifee_pressure_constr2 = Constraint(expr=
                m.fs.menifee_a.feed.pressure[0] == m.fs.menifee_b.feed.pressure[0])
            m.fs.perris_pressure_constr1 = Constraint(expr=
                m.fs.perris_i_a.feed.pressure[0] <= 14)
            m.fs.perris_pressure_constr2 = Constraint(expr=
                m.fs.perris_i_a.feed.pressure[0] == m.fs.perris_i_b.feed.pressure[0])
            m.fs.perris_recov_constr1 = Constraint(expr=
                m.fs.perris_i_a.flow_vol_out[0] / m.fs.perris_i_a.flow_vol_in[0] <= 0.75)
            m.fs.perris_recov_constr1 = Constraint(expr=
                m.fs.perris_i_a.flow_vol_out[0] / m.fs.perris_i_a.flow_vol_in[0] >= 0.70)
            m.fs.area_constr1 = Constraint(expr=
                (m.fs.perris_i_a.membrane_area[0] + m.fs.perris_i_b.membrane_area[0]) / 
                (m.fs.menifee_a.membrane_area[0] + m.fs.menifee_b.membrane_area[0]) >= 1.4)
            m.fs.area_constr2 = Constraint(expr=
                (m.fs.perris_i_a.membrane_area[0] + m.fs.perris_i_b.membrane_area[0]) / 
                (m.fs.menifee_a.membrane_area[0] + m.fs.menifee_b.membrane_area[0]) <= 1.6)

        elif 'zld' in scenario:
            m.fs.first_pass_press_constr1 = Constraint(expr=
                m.fs.menifee_first_pass.feed.pressure[0] <= 14)
            m.fs.first_pass_press_constr2 = Constraint(expr=
                m.fs.menifee_first_pass.feed.pressure[0] == m.fs.perris_i_first_pass.feed.pressure[0])
            m.fs.second_pass_press_constr1 = Constraint(expr=
                m.fs.menifee_second_pass.feed.pressure[0] <= 18)
            m.fs.second_pass_press_constr2 = Constraint(expr=
                m.fs.menifee_second_pass.feed.pressure[0] == m.fs.perris_i_second_pass.feed.pressure[0])
            m.fs.area_constr1 = Constraint(expr=
                (m.fs.perris_i_first_pass.membrane_area[0] + 
                m.fs.perris_i_second_pass.membrane_area[0]) / 
                (m.fs.menifee_first_pass.membrane_area[0] + 
                m.fs.menifee_second_pass.membrane_area[0]) >= 1.4)
            m.fs.area_constr2 = Constraint(expr=
                (m.fs.perris_i_first_pass.membrane_area[0] + 
                m.fs.perris_i_second_pass.membrane_area[0]) / 
                (m.fs.menifee_first_pass.membrane_area[0] + 
                m.fs.menifee_second_pass.membrane_area[0]) <= 1.6)
            m.fs.ro_area_constr1 = Constraint(expr=
                m.fs.menifee_first_pass.membrane_area[0] <= m.fs.menifee_second_pass.membrane_area[0])
            m.fs.ro_area_constr2 = Constraint(expr=
                m.fs.perris_i_first_pass.membrane_area[0] <= m.fs.perris_i_second_pass.membrane_area[0])
            m.fs.area_ratio_constr1 = Constraint(expr=
                (m.fs.menifee_first_pass.membrane_area[0] / 
                m.fs.menifee_second_pass.membrane_area[0]) == 
                (m.fs.perris_i_first_pass.membrane_area[0] /
                 m.fs.perris_i_second_pass.membrane_area[0]))

    if case_study == 'ocwd':  # Facility data in email from Dan Giammar 7/7/2021
        m.fs.ro_pressure_constr = Constraint(expr=
            m.fs.reverse_osmosis.feed.pressure[0] <= 15)  # Facility data: RO pressure is 140-220 psi (~9.7-15.1 bar)
        m.fs.microfiltration.water_recovery.fix(0.9)

    if case_study == 'produced_water_injection' and scenario == 'swd_well':
        m.fs.brine_concentrator.water_recovery.fix(0.725)

    if case_study == "uranium":
        m.fs.ion_exchange.removal_fraction[0, "tds"].unfix()
        m.fs.ion_exchange.water_recovery.fix(0.967)
        m.fs.ion_exchange.anion_res_capacity.unfix()
        m.fs.ion_exchange.cation_res_capacity.unfix()

    if case_study == 'irwin':
        m.fs.brine_concentrator.water_recovery.fix(0.8)

    return m