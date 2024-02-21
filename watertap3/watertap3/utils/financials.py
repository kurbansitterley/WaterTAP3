##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
##############################################################################

import os
import pandas as pd
from pyomo.environ import (Block, Expression, Constraint, Param, Var, NonNegativeReals, units as pyunits)

from .ml_regression import get_linear_regression

__all__ = ['SystemSpecs', 'get_complete_costing', 'get_ind_table', 'get_system_specs',
           'get_system_costing', 'global_costing_parameters']

last_year_for_cost_indicies = 2050

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
case_study_TEA_basis_file = os.path.abspath(os.path.join(__location__, os.pardir)) + "/data/case_study_TEA_basis.csv"
industrial_electricity_costs_file = os.path.abspath(os.path.join(__location__, os.pardir)) + "/data/industrial_electricity_costs_2020.csv"
chem_costs_file = os.path.abspath(os.path.join(__location__, os.pardir)) + "/data/chemical_costs.csv"
plant_costs_indices_file = os.path.abspath(os.path.join(__location__, os.pardir)) + "/data/plant_cost_indices.csv"

class SystemSpecs():

    def __init__(self, train=None):
        basis_data = pd.read_csv(case_study_TEA_basis_file, index_col='case_study')
        elec_cost = pd.read_csv(industrial_electricity_costs_file, index_col='location')
        elec_cost.index = elec_cost.index.str.lower()
        case_study = train['case_study']
        if 'test' in case_study:
            case_study = 'test'
        self.location = basis_data[basis_data['variable'] == 
            'location_basis'].loc[case_study].value
        self.elec_price = float(elec_cost.loc[self.location])
        self.land_cost_percent_FCI = float(basis_data[basis_data['variable'] == 
            'land_cost_percent'].loc[case_study].value)
        self.working_cap_percent_FCI = float(basis_data[basis_data['variable'] == 
            'working_capital_percent'].loc[case_study].value)
        self.salaries_percent_FCI = float(basis_data[basis_data['variable'] == 
            'base_salary_per_fci'].loc[case_study].value)
        self.maintenance_costs_percent_FCI = float(basis_data[basis_data['variable'] == 
            'maintenance_cost_percent'].loc[case_study].value)
        self.lab_fees_percent_FCI = float(basis_data[basis_data['variable'] == 
            'laboratory_fees_percent'].loc[case_study].value)
        self.insurance_taxes_percent_FCI = float(basis_data[basis_data['variable'] == 
            'insurance_and_taxes_percent'].loc[case_study].value)
        self.benefit_percent_of_salary = float(basis_data[basis_data['variable'] == 
            'employee_benefits_percent'].loc[case_study].value)
        self.plant_lifetime_yrs = int(basis_data[basis_data['variable'] == 
            'plant_life_yrs'].loc[case_study].value)
        self.analysis_yr_cost_indices = int(basis_data[basis_data['variable'] == 
            'analysis_year'].loc[case_study].value)
        self.debt_interest_rate = float(basis_data[basis_data['variable'] == 
            'debt_interest_rate'].loc[case_study].value)
        self.plant_cap_utilization = float(basis_data[basis_data['variable'] == 
            'plant_cap_utilization'].loc[case_study].value)


def get_complete_costing(costing, basis_year=2020, tpec_tic=None):
    '''
    Function to build costing block for each WaterTAP3 unit.

    :param costing: Costing block object from WaterTAP3 unit model.
    :type costing: object
    :return:
    '''
    unit = costing.parent_block()
    costing.basis_year = basis_year
    sys_cost_params = unit.parent_block().costing_param
    time = unit.flowsheet().config.time
    t = time.first()
    flow_in_m3yr = pyunits.convert(unit.flow_vol_in[t], 
            to_units=pyunits.m**3/pyunits.year)

    if tpec_tic == 'TPEC':
        unit.tpec_tic.fix(sys_cost_params.tpec)
    elif tpec_tic == 'TIC':
        unit.tpec_tic.fix(sys_cost_params.tic)

    basis_year = costing.basis_year
    sys_specs = unit.parent_block().costing_param

    chem_dict = unit.chem_dict

    ## COSTING INDICES
    df = get_ind_table(sys_specs.analysis_yr_cost_indices)
    costing.ind_df = df.copy()
    costing.cap_replacement_parts = df.loc[basis_year].Capital_Factor
    costing.catalysts_chemicals = df.loc[basis_year].CatChem_Factor
    costing.labor_and_other_fixed = df.loc[basis_year].Labor_Factor
    costing.consumer_price_index = df.loc[basis_year].CPI_Factor

    costing.fixed_cap_inv = (costing.fixed_cap_inv_unadjusted * costing.cap_replacement_parts) 
    costing.land_cost = costing.fixed_cap_inv * sys_specs.land_cost_percent_FCI
    costing.working_cap = costing.fixed_cap_inv * sys_specs.working_cap_percent_FCI
    costing.contingency = costing.fixed_cap_inv * sys_specs.contingency_cost_percent_FCI
    costing.component_replacement = costing.fixed_cap_inv * sys_specs.component_replace_percent_FCI
    costing.base_employee_salary_cost = costing.fixed_cap_inv_unadjusted * sys_specs.salaries_percent_FCI
    costing.salaries = costing.labor_and_other_fixed * costing.base_employee_salary_cost
    costing.benefits = costing.salaries * sys_specs.benefit_percent_of_salary
    costing.maintenance = costing.fixed_cap_inv * sys_specs.maintenance_costs_percent_FCI
    costing.lab = costing.fixed_cap_inv * sys_specs.lab_fees_percent_FCI
    costing.insurance_taxes = costing.fixed_cap_inv * sys_specs.insurance_taxes_percent_FCI

    cat_chem_df = pd.read_csv(chem_costs_file, index_col='Material')
    costing.chem_cost_sum = chem_cost_sum = 0
    for chem, dose in chem_dict.items():
        # if chem == 'unit_cost':
        #     chem_cost_sum = dose * costing.fixed_cap_inv * 1E6
        # else:
        chem_name = chem.replace('(', '').replace(')', '').replace(' ', '_').replace('%', 'pct')
        setattr(costing, f'{chem_name}_unit_price', \
            Var(initialize=0.1,
                bounds=(0, None),
                doc=f'Unit Cost of {chem}'))
        setattr(costing, f'{chem_name}_dose', \
            Var(initialize=0.1,
                bounds=(0, None),
                doc=f'Dose of {chem}'))
        chem_price_var = getattr(costing, f'{chem_name}_unit_price')
        chem_dose_var = getattr(costing, f'{chem_name}_dose')
        
        if unit.unit_type == 'chlorination':
            chem_dose_var.fix(dose)
        elif chem == 'unit_cost':
            chem_dose_var.fix(dose)
            chem_cost_sum += dose * costing.fixed_cap_inv * 1E6
        else:
            chem_price_var.fix(cat_chem_df.loc[chem].Price)
            chem_dose_var.fix(dose())
        # chem_cost_sum += costing.cap_replacement_parts * flow_in_m3yr * \
        #                 chem_price_var * dose/0.3 * sys_specs.plant_cap_utilization
        chem_cost_sum += costing.catalysts_chemicals * flow_in_m3yr * \
                        chem_price_var * dose * sys_specs.plant_cap_utilization
    # costing.chem_cost_sum = chem_cost_sum
    costing.cat_and_chem_cost = Var()
    
    costing.cat_and_chem_cost_costr = Constraint(expr=costing.cat_and_chem_cost == ((chem_cost_sum * 1E-6) * 
            (1 - costing.catchem_reduction)) * costing.catchem_uncertainty)

    # if not hasattr(costing, 'electricity_cost'):
    costing.electricity_intensity = (unit.electricity * 
        (1 - costing.elect_intens_reduction)) * costing.elect_intens_uncertainty
    costing.electricity_cost = ((costing.electricity_intensity * \
        flow_in_m3yr * sys_specs.electricity_price * 1E-6) * 
            sys_specs.plant_cap_utilization) 

    if not hasattr(costing, 'other_var_cost'):
        costing.other_var_cost = 0

    else:
        costing.other_var_cost = (costing.other_var_cost * \
            (1 - costing.other_reduction)) * costing.other_uncertainty

    costing.total_cap_investment = (costing.fixed_cap_inv + costing.land_cost + costing.working_cap)
    costing.total_fixed_op_cost = ((costing.salaries + costing.benefits + 
            costing.maintenance + costing.lab + costing.insurance_taxes))
    costing.annual_op_main_cost = ((costing.cat_and_chem_cost + costing.electricity_cost + 
            costing.other_var_cost + costing.total_fixed_op_cost)) 
    costing.total_operating_cost = ((costing.total_fixed_op_cost + costing.cat_and_chem_cost + 
            costing.electricity_cost + costing.other_var_cost) )


def get_ind_table(analysis_yr_cost_indices):
    '''
    Function to get costing indicies for WaterTAP3 model.

    :param analysis_yr_cost_indices: Year to get costing indices for.
    :type analysis_yr_cost_indices: int
    :return: Indicies DataFrame
    '''
    df = pd.read_csv(plant_costs_indices_file)

    df1 = pd.DataFrame()
    for name in df.columns[1:]:
        a, b = get_linear_regression(list(df.year), list(df[(f'{name}')]), name)
        new_list = []
        yr_list = []
        for yr in range(df.year.max() + 1, last_year_for_cost_indicies + 1):
            new_list.append(a * yr + b)
            yr_list.append(yr)
        df1[name] = new_list
    df1['year'] = yr_list
    df = pd.concat([df, df1], axis=0)

    new_cost_variables = ['capital', 'chemical', 'labor']
    for variable in new_cost_variables:
        ind_name = '%s_index' % variable
        fac_name = '%s_factor' % variable
        df[fac_name] = (df[df.year == analysis_yr_cost_indices][ind_name].max() / df[ind_name])
    df = df.set_index(df.year)
    df = df.replace(1.0, 1.00000000001)

    return df


def get_system_specs(m_fs):
    '''
    Function to set costing parameters for WaterTAP3 model.

    '''
    m_fs.costing_param = Block()
    b = m_fs.costing_param

    b.electricity_price = Var(
        initialize=0.07,
        doc='Electricity cost [$/kWh]')
    b.maintenance_costs_percent_FCI = Var(
        initialize=0.008,
        doc='Maintenance/contingency cost as % FCI')
    b.salaries_percent_FCI = Var(
        initialize=0.001,
        doc='Salaries cost as % FCI')
    b.benefit_percent_of_salary = Var(
        initialize=0.07,
        doc='Benefits cost as % FCI')
    b.insurance_taxes_percent_FCI = Var(
        initialize=0.002,
        doc='Insurance/taxes cost as % FCI')
    b.lab_fees_percent_FCI = Var(
        initialize=0.003,
        doc='Lab cost as % FCI')
    b.land_cost_percent_FCI = Var(
        initialize=0.0015,
        doc='Land cost as % FCI')
    b.plant_lifetime_yrs = Var(
        initialize=30,
        doc='Plant lifetime [years]')
    b.plant_cap_utilization = Var(
        initialize=1,
        doc='Plant capacity utilization [%]')
    b.working_cap_percent_FCI = Var(
        initialize=0.008,
        doc='Working capital as % FCI')
    b.wacc = Var(
        initialize=0.05,
        doc='Weighted Average Cost of Capital (WACC)')
    b.contingency_cost_percent_FCI = Var(
        initialize=0,
        doc='Contingency costs as % FCI')
    b.component_replace_percent_FCI = Var(
        initialize=0,
        doc='Component replacement costs as % FCI')

    system_specs = SystemSpecs(m_fs.train)

    b.electricity_price.fix(system_specs.elec_price)
    b.salaries_percent_FCI.fix(system_specs.salaries_percent_FCI)
    b.land_cost_percent_FCI.fix(system_specs.land_cost_percent_FCI)
    b.maintenance_costs_percent_FCI.fix(system_specs.maintenance_costs_percent_FCI)
    b.lab_fees_percent_FCI.fix(system_specs.lab_fees_percent_FCI)
    b.insurance_taxes_percent_FCI.fix(system_specs.insurance_taxes_percent_FCI)
    b.plant_lifetime_yrs.fix(system_specs.plant_lifetime_yrs)

    b.benefit_percent_of_salary.fix(system_specs.benefit_percent_of_salary)
    b.working_cap_percent_FCI.fix(system_specs.working_cap_percent_FCI)
    b.plant_cap_utilization.fix(system_specs.plant_cap_utilization)  # 1.0
    b.wacc.fix(system_specs.debt_interest_rate)
    b.contingency_cost_percent_FCI.fix(0)
    b.component_replace_percent_FCI.fix(0)

    b.analysis_yr_cost_indices = system_specs.analysis_yr_cost_indices
    b.location = system_specs.location

    b.tpec = 3.4
    b.tic = 1.65


def get_system_costing(m_fs):
    '''
    Function to aggregate unit model results for calculation of system costing for WaterTAP3 model.

    '''
    if not hasattr(m_fs, 'costing'):
        m_fs.costing = Block()
    b = m_fs.costing
    time = m_fs.config.time
    t = time.first()
    sys_specs = m_fs.costing_param

    total_capital_investment_var_lst = []
    cat_and_chem_cost_lst = []
    electricity_cost_lst = []
    other_var_cost_lst = []
    total_fixed_op_cost_lst = []
    electricity_intensity_lst = []

    wacc = sys_specs.wacc

    b.capital_recovery_factor = (wacc * (1 + wacc) ** sys_specs.plant_lifetime_yrs) / (
            ((1 + wacc) ** sys_specs.plant_lifetime_yrs) - 1)

    for b_unit in m_fs.component_objects(Block, descend_into=True):
        if hasattr(b_unit, 'costing'):
            total_capital_investment_var_lst.append(b_unit.costing.total_cap_investment)
            cat_and_chem_cost_lst.append(b_unit.costing.cat_and_chem_cost)
            electricity_cost_lst.append(b_unit.costing.electricity_cost)
            other_var_cost_lst.append(b_unit.costing.other_var_cost)
            total_fixed_op_cost_lst.append(b_unit.costing.total_fixed_op_cost)

    b.cat_and_chem_cost_annual = Expression(expr=
        (sum(cat_and_chem_cost_lst) * (1 - b.sys_catchem_reduction)) * b.sys_catchem_uncertainty)
    b.electricity_cost_annual = Expression(expr=
        (sum(electricity_cost_lst) * (1 - b.sys_elect_reduction)) * b.sys_elect_uncertainty)
    b.other_var_cost_annual = Expression(expr=
        (sum(other_var_cost_lst) * (1 - b.sys_other_reduction)) * b.sys_other_uncertainty)
    b.fixed_op_cost_annual = Expression(expr=
        (sum(total_fixed_op_cost_lst) * (1 - b.sys_fixed_op_reduction)) * b.sys_fixed_op_uncertainty)
    b.operating_cost_annual = Expression(expr=
        (b.fixed_op_cost_annual + b.cat_and_chem_cost_annual + 
        b.electricity_cost_annual + b.other_var_cost_annual))
    #
    b.capital_investment_total = Expression(expr=
        (sum(total_capital_investment_var_lst) * (1 - b.sys_tci_reduction)) * b.sys_tci_uncertainty)
    b.cat_and_chem_cost_total = Expression(expr=
        b.cat_and_chem_cost_annual * m_fs.costing_param.plant_lifetime_yrs)
    b.electricity_cost_total = Expression(expr=
        b.electricity_cost_annual * m_fs.costing_param.plant_lifetime_yrs)
    b.other_var_cost_total = Expression(expr=
        b.other_var_cost_annual * m_fs.costing_param.plant_lifetime_yrs)
    b.fixed_op_cost_total = Expression(expr=
        b.fixed_op_cost_annual * m_fs.costing_param.plant_lifetime_yrs)
    b.operating_cost_total = Expression(expr=
        ((b.fixed_op_cost_total + b.cat_and_chem_cost_total + b.electricity_cost_total + 
        b.other_var_cost_total) * (1 - b.sys_total_op_reduction)) * b.sys_total_op_uncertainty)

    recovered_water_flow = 0

    time = m_fs.config.time.first()

    for b_unit in m_fs.component_objects(Block, descend_into=False):
        if hasattr(b_unit, 'outlet'):
            if len(getattr(b_unit, 'outlet').arcs()) == 0:
                if hasattr(b_unit.parent_block(), 'pfd_dict'):
                    if b_unit.parent_block().pfd_dict[str(b_unit)[3:]]['Type'] == 'use':
                        recovered_water_flow += b_unit.flow_vol_out[time]
                else:
                    if 'reverse_osmosis' in str(b_unit):
                        recovered_water_flow += b_unit.flow_vol_out[time]
                    if 'cooling_tower' in str(b_unit):
                        recovered_water_flow += b_unit.flow_vol_out[time]

    b.treated_water = recovered_water_flow
    intake_str = [k for k, v in m_fs.pfd_dict.items() if v['Type'] == 'intake']
    waste_str = [k for k, v in m_fs.pfd_dict.items() if v['Type'] == 'waste']

    intakes = [unit for unit in m_fs.component_objects(Block, descend_into=False) if
               hasattr(unit, 'flow_vol_in') and str(unit)[3:] in intake_str]
    wastes = [unit for unit in m_fs.component_objects(Block, descend_into=False) if
              hasattr(unit, 'flow_vol_in') and str(unit)[3:] in waste_str]

    flow_ins = [getattr(intake, 'flow_vol_in')[0]() for intake in intakes]
    flow_wastes = [getattr(waste, 'flow_vol_in')[0]() for waste in wastes]

    m_fs.sys_flow_in = sum(flow_ins)
    m_fs.sys_flow_waste = sum(flow_wastes)

    b.sum_of_inflow = sum_of_inflow = 0
    for key in b.parent_block().flow_in_dict.keys():
        sum_of_inflow += getattr(m_fs, key).flow_vol_in[time]

    b.system_recovery = b.treated_water / sum_of_inflow

    # LCOW for each unit
    
    for b_unit in m_fs.component_objects(Block, descend_into=True):
        if hasattr(b_unit, 'costing'):
            setattr(b_unit, 'LCOW', Expression(
                    expr=1E6 * (b_unit.costing.total_cap_investment * b.capital_recovery_factor + b_unit.costing.annual_op_main_cost) /
                         (b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization),
                    doc='Unit Levelized Cost of Water [$/m3]'))

            setattr(b_unit, 'LCOW_TCI', Expression(
                    expr=1E6 * (b_unit.costing.total_cap_investment * b.capital_recovery_factor) /
                         (b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization),
                    doc='Unit TCI Levelized Cost of Water [$/m3]'))

            setattr(b_unit, 'LCOW_elec', Expression(
                    expr=1E6 * (b_unit.costing.electricity_cost) /
                         (b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization),
                    doc='Unit Electricity Levelized Cost of Water [$/m3]'))

            setattr(b_unit, 'LCOW_fixed_op', Expression(
                    expr=1E6 * (b_unit.costing.total_fixed_op_cost) /
                         (b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization),
                    doc='Unit Fixed Operating Levelized Cost of Water [$/m3]'))

            setattr(b_unit, 'LCOW_chem', Expression(
                    expr=1E6 * (b_unit.costing.cat_and_chem_cost) /
                         (b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization),
                    doc='Unit Chemical Levelized Cost of Water [$/m3]'))

            setattr(b_unit, 'LCOW_other', Expression(
                    expr=1E6 * (b_unit.costing.other_var_cost) /
                         (b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization),
                    doc='Unit Other O&M Levelized Cost of Water [$/m3]'))

            setattr(b_unit, 'LCOW_total_op', Expression(
                    expr=1E6 * (b_unit.costing.total_operating_cost) /
                         (b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization),
                    doc='Unit Total Operating Levelized Cost of Water [$/m3]'))

            setattr(b_unit, 'elec_int_treated', Expression(
                    expr=(b_unit.costing.electricity_cost * 1E6 / sys_specs.electricity_price) /
                         (b.treated_water * 3600 * 24 * 365),
                    doc='Unit Electricity Intensity [kWh/m3]'))

    b.LCOW_TCI = Expression(expr=1E6 * (b.capital_investment_total * b.capital_recovery_factor) / (
            b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization))

    b.LCOW_elec = Expression(expr=1E6 * (b.electricity_cost_annual) / (
            b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization))

    b.LCOW_fixed_op = Expression(expr=1E6 * (b.fixed_op_cost_annual) / (
            b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization))

    b.LCOW_chem = Expression(expr=1E6 * (b.cat_and_chem_cost_annual) / (
            b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization))

    b.LCOW_other_onm = Expression(expr=1E6 * (b.other_var_cost_annual) / (
            b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization))

    b.LCOW_total_op = Expression(expr=1E6 * (b.operating_cost_annual) / (
            b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization))

    ## GET TOTAL ELECTRICITY CONSUMPTION IN kWh/m3 of treated water
    b.electricity_intensity = Expression(
            expr=(b.electricity_cost_annual * 1E6 / sys_specs.electricity_price) /
                 (b.treated_water * 3600 * 24 * 365),
            doc='Electricity Intensity [kWh/m3]')

    b.LCOW = Expression(
            expr=1E6 * (b.capital_investment_total * b.capital_recovery_factor + b.operating_cost_annual) /
                 (b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization),
            doc='Levelized Cost of Water [$/m3]')

    b.LCOW_inflow = Expression(
            expr=1E6 * (b.capital_investment_total * b.capital_recovery_factor + b.operating_cost_annual) /
                 (sum_of_inflow * 3600 * 24 * 365 * sys_specs.plant_cap_utilization),
            doc='Levelized Cost of Water by influent flow [$/m3]')

    b.elec_frac_LCOW = Expression(
            expr=((1E6 * (b.electricity_cost_annual) /
                   (b.treated_water * 3600 * 24 * 365 * sys_specs.plant_cap_utilization))) / b.LCOW,
            doc='Electricity cost as fraction of LCOW')

def global_costing_parameters(self, year=None):
    if year is None:
        year = '2018'
    ce_index_dic = {
            '2019': 680,
            '2018': 671.1,
            '2017': 567.5,
            '2016': 541.7,
            '2015': 556.8,
            '2014': 576.1,
            '2013': 567.3,
            '2012': 584.6,
            '2011': 585.7,
            '2010': 550.8
            }
