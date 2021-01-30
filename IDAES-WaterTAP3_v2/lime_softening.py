##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Demonstration zeroth-order model for WaterTAP3
"""

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import Block, Constraint, Var, units as pyunits
from pyomo.network import Port

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block

from pyomo.environ import (
    Expression, Var, Param, NonNegativeReals, units as pyunits)

# Import WaterTAP# financials module
import financials
from financials import *

from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory
from pyomo.network import Arc
from idaes.core import FlowsheetBlock

# Import properties and units from "WaterTAP Library"
from water_props import WaterParameterBlock

# Set inlet conditions to first unit
# IDAES Does have Feed components for this, but that would require a bit
# more work to set up to work (as it relies on things in the property package
# that aren't implemented for this example).
# I am just picking numbers for most of these


### FACTORS FOR ZEROTH ORDER MODEL -> TODO -> READ IN AUTOMATICALLY BASED ON UNIT PROCESS --> CREATE TABLE?!###
flow_recovery_factor = 0.99999
tds_removal_factor = 0

# Perfomance Parameter Values for Process: Constituent removals.
toc_removal_factor = 0.0  
nitrates_removal_factor = 0.0  
TOrC_removal = 0.0  
EEQ_removal = 0.0 
NDMA_removal = 0.0
PFOS_PFOA_removal = 0.0
protozoa_removal = 0.0
virus_removal = 0.0

# capital costs basis

base_fixed_cap_cost = 2.8403 # from VAR tab
cap_scaling_exp = .7306  # from VAR tab

basis_year = 2020
fixed_op_cost_scaling_exp = 0.7


# recycle_factor = (1 - recovery_factor) * (recyle_fraction_of_waste)
waste_factor = 1 - flow_recovery_factor  # - G.edges[edge]['recycle_factor']


### data specific to UP - should probably be a chemical addition class ###
### TODO check these against original source: 
#Minnesota Rural Water Association, Chapter 16 Lime Softening #(https://www.mrwa.com/WaterWorksMnl/Chapter%2016%20Lime%20Softening.pdf)
#https://www.necoindustrialwater.com/analysis-ion-exchange-vs-lime-softening/

#### CAT AND CHEMS START ###
#these are the constituents in the inlet water assuming seawater 35 g/L -> ASSUMED TO ONLY BE VALID FOR CARLSBAD?
co2 = 1.27 # mg/L NOT USED
bicarb_alk_CaCO3 = 0.56 # mg/L NOT USED
hydrox_alk_CaCO3 = 0.56 # mg/L NOT USED
magnesium_dissolved_lime = 2.3 # mg/L
magnesium_dissolved_factor = 30 # TODO
lime = magnesium_dissolved_factor * magnesium_dissolved_lime / 1000 
# assumed to be Lime_Suspension_(Ca(OH)2). kg/m3
Lime_Suspension_CaOH_2 = lime / 264.172 # converted to kg/gal todo in pyunits

chem_dic = {"Lime_Suspension_CaOH_2" : Lime_Suspension_CaOH_2}

#### CAT AND CHEMS END ###

# You don't really want to know what this decorator does
# Suffice to say it automates a lot of Pyomo boilerplate for you
@declare_process_block_class("UnitProcess")
class UnitProcessData(UnitModelBlockData):
       
    """
    This class describes the rules for a zeroth-order model for a unit
    """
    # The Config Block is used tpo process arguments from when the model is
    # instantiated. In IDAES, this serves two purposes:
    #     1. Allows us to separate physical properties from unit models
    #     2. Lets us give users options for configuring complex units
    # For WaterTAP3, this will mainly be boilerplate to keep things consistent
    # with ProteusLib and IDAES.
    # The dynamic and has_holdup options are expected arguments which must exist
    # The property package arguments let us define different sets of contaminants
    # without needing to write a new model.
    CONFIG = ConfigBlock()
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not,
**default** = False. Equilibrium Reactors do not support dynamic behavior."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Indicates whether holdup terms should be constructed or not.
**default** - False. Equilibrium reactors do not have defined volume, thus
this must be False."""))
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))
    
    from unit_process_equations import initialization
    #unit_process_equations.get_base_unit_process()

    #build(up_name = "nanofiltration_twb")
    
    def build(self):
        import unit_process_equations
        return unit_process_equations.build_up(self, up_name_test = "lime_softening")
    
    
    def get_costing(self, module=financials, cost_method="wt", year=None):
        """
        We need a get_costing method here to provide a point to call the
        costing methods, but we call out to an external consting module
        for the actual calculations. This lets us easily swap in different
        methods if needed.

        Within IDAES, the year argument is used to set the initial value for
        the cost index when we build the model.
        """
        # First, check to see if global costing module is in place
        # Construct it if not present and pass year argument
        if not hasattr(self.flowsheet(), "costing"):
            self.flowsheet().get_costing(module=module, year=year)

        # Next, add a sub-Block to the unit model to hold the cost calculations
        # This is to let us separate costs from model equations when solving
        self.costing = Block()
        # Then call the appropriate costing function out of the costing module
        # The first argument is the Block in which to build the equations
        # Can pass additional arguments as needed
        
        #up_costing(self.costing, cost_method=cost_method)
        
        # There are a couple of variables that IDAES expects to be present
        # These are fairly obvious, but have pre-defined names
        def _make_vars(self):
            # build generic costing variables (all costing models need these vars)
            self.base_cost = Var(initialize=1e5,
                                 domain=NonNegativeReals,
                                 doc='Unit Base Cost cost in $')
            self.purchase_cost = Var(initialize=1e4,
                                     domain=NonNegativeReals,
                                     doc='Unit Purchase Cost in $')
    
    
        # Build a costing method for each type of unit
        def up_costing(self, cost_method="wt"):
            
            '''
            This is where you create the variables and equations specific to each unit.
            This method should mainly consider capital costs for the unit - operating
            most costs should done for the entire flowsheet (e.g. common utilities).
            Unit specific operating costs, such as chemicals, should be done here with
            standard names that can be collected at the flowsheet level.

            You can access variables from the unit model using:

                self.parent_block().variable_name

            You can also have unit specific parameters here, which could be retrieved
            from the spreadsheet
            '''
            
            # Lime softening unit process based on Texas Water Development Board IT3PR.
            # User Manual ==> https://www.twdb.texas.gov/publications/reports/contracted_reports/doc/1348321632_manual.pdf
            
            
    
            _make_vars(self)

            self.base_fixed_cap_cost = Param(mutable=True,
                                             initialize=base_fixed_cap_cost,
                                             doc="Some parameter from TWB")
            self.cap_scaling_exp = Param(mutable=True,
                                         initialize=cap_scaling_exp,
                                         doc="Another parameter from TWB")

            # Get the first time point in the time domain
            # In many cases this will be the only point (steady-state), but lets be
            # safe and use a general approach
            time = self.parent_block().flowsheet().config.time.first()

            # Get the inlet flow to the unit and convert to the correct units
            flow_in = pyunits.convert(self.parent_block().flow_vol_in[time],
                                      to_units=pyunits.Mgallons/pyunits.day)
            

            ################### TWB METHOD ###########################################################
            if cost_method == "twb":
                    self.fixed_cap_inv_unadjusted = Expression(
                        expr=self.base_fixed_cap_cost *
                        flow_in ** self.cap_scaling_exp,
                        doc="Unadjusted fixed capital investment")
            ##############################################################################

            ################## WATERTAP METHOD ###########################################################
            if cost_method == "wt":

                # cost index values - TODO MOVE THIS TO TOP
                df = get_ind_table()
                self.cap_replacement_parts = df.loc[basis_year].Capital_Factor
                self.catalysts_chemicals = df.loc[basis_year].CatChem_Factor
                self.labor_and_other_fixed = df.loc[basis_year].Labor_Factor
                self.consumer_price_index = df.loc[basis_year].CPI_Factor

                # capital costs (unit: MM$) ---> TCI IN EXCEL
                self.fixed_cap_inv_unadjusted = Expression(
                    expr=self.base_fixed_cap_cost * flow_in ** self.cap_scaling_exp,
                    doc="Unadjusted fixed capital investment") 

                self.fixed_cap_inv = self.fixed_cap_inv_unadjusted * self.cap_replacement_parts
                self.land_cost = self.fixed_cap_inv * land_cost_precent_FCI
                self.working_cap = self.fixed_cap_inv * working_cap_precent_FCI
                self.total_cap_investment = self.fixed_cap_inv + self.land_cost + self.working_cap

                # variable operating costs (unit: MM$/yr) -> MIKE TO DO -> ---> CAT+CHEM IN EXCEL
                # --> should be functions of what is needed!?
                # cat_chem_df = pd.read_csv('catalyst_chemicals.csv')
                # cat_and_chem = flow_in * 365 * on_stream_factor # TODO
                self.electricity = .01  # kwh/m3 given in PML tab, no source TODO
                
                # TODO for systematic approach elsewhere
                cat_chem_df = pd.read_csv('data/catalyst_chemicals.csv', index_col = "Material")
                chem_cost_sum = 0 
                for key in chem_dic.keys():
                    chem_cost = cat_chem_df.loc[key].Price
                    chem_cost_sum = chem_cost_sum + self.catalysts_chemicals * flow_in * 365 * chem_cost * chem_dic[key] * on_stream_factor
                self.cat_and_chem_cost = chem_cost_sum
                
                
                flow_in_m3yr = (pyunits.convert(self.parent_block().flow_vol_in[time],
                                      to_units=pyunits.m**3/pyunits.year))
                self.electricity_cost = Expression(
                        expr= (self.electricity * flow_in_m3yr * elec_price/1000000),
                        doc="Electricity cost") # M$/yr
                
                self.other_var_cost = 0 # Expression(
                        #expr= self.cat_and_chem_cost - self.electricity_cost,
                        #doc="Other variable cost")

                # fixed operating cost (unit: MM$/yr)  ---> FIXED IN EXCEL
#                 self.base_employee_salary_cost = self.base_fixed_cap_cost * flow_in ** self.cap_scaling_exp * salaries_percent_FCI
#                 self.salaries = (
#                     self.labor_and_other_fixed
#                     * self.base_employee_salary_cost
#                     * flow_in ** fixed_op_cost_scaling_exp
#                 )
                
#                 self.salaries = (
#                     (self.labor_and_other_fixed ** fixed_op_cost_scaling_exp) * (salaries_percent_FCI 
#                           * self.fixed_cap_inv_unadjusted) ** fixed_op_cost_scaling_exp)
               
                self.base_employee_salary_cost = self.fixed_cap_inv_unadjusted * salaries_percent_FCI
                self.salaries = Expression(
                        expr= self.labor_and_other_fixed * self.base_employee_salary_cost,
                        doc="Salaries")
                
                self.benefits = self.salaries * benefit_percent_of_salary
                self.maintenance = maintinance_costs_precent_FCI * self.fixed_cap_inv
                self.lab = lab_fees_precent_FCI * self.fixed_cap_inv
                self.insurance_taxes = insurance_taxes_precent_FCI * self.fixed_cap_inv
                self.total_fixed_op_cost = Expression(
                    expr = self.salaries + self.benefits + self.maintenance + self.lab + self.insurance_taxes)

                self.total_up_cost = (
                    self.total_cap_investment
                    + self.cat_and_chem_cost
                    + self.electricity_cost
                    + self.other_var_cost
                    + self.total_fixed_op_cost
                
                )

            #return total_up_cost
    
        up_costing(self.costing, cost_method=cost_method)
          
        
# OTHER CALCS

def create(m, up_name):
    
    # Set removal and recovery fractions
    getattr(m.fs, up_name).water_recovery.fix(flow_recovery_factor)
    getattr(m.fs, up_name).removal_fraction[:, "TDS"].fix(tds_removal_factor)
    # I took these values from the WaterTAP3 nf model
    getattr(m.fs, up_name).removal_fraction[:, "TOC"].fix(toc_removal_factor)
    getattr(m.fs, up_name).removal_fraction[:, "nitrates"].fix(nitrates_removal_factor)

    # Also set pressure drops - for now I will set these to zero
    getattr(m.fs, up_name).deltaP_outlet.fix(1e-4)
    getattr(m.fs, up_name).deltaP_waste.fix(1e-4)

    # Adding costing for units - this is very basic for now so use default settings
    getattr(m.fs, up_name).get_costing(module=financials)

    return m        
        
        
           
        
        
        