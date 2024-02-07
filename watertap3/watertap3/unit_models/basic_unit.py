from pyomo.environ import Expression, units as pyunits
from watertap3.utils import cost_curves, financials
from watertap3.core.wt3_unit import WT3UnitProcess

## REFERENCE: ADD REFERENCE HERE

module_name = 'basic_unit'

class UnitProcess(WT3UnitProcess):

    def fixed_cap(self):
        '''
        :param flow_in: Flow in to basic unit [m3/hr]
        :type flow_in: float
        '''
        time = self.flowsheet().config.time.first()
        sys_cost_params = self.parent_block().costing_param
        flow_in_m3yr = pyunits.convert(self.flow_in,
            to_units=(pyunits.m**3 / pyunits.year))

        if self.unit_process_name == "tramp_oil_tank":
            disposal_cost = 0.000114 # Kiran's disposal cost assumption $/m3
            self.costing.other_var_cost = flow_in_m3yr * disposal_cost * sys_cost_params.plant_cap_utilization * 1E-6

        if self.kind == 'flow':
            flow_basis = self.basis * (pyunits.m**3/pyunits.hour)
            flow_factor = self.flow_in / flow_basis
            basic_cap = self.cap_basis * flow_factor ** self.cap_exp
            return basic_cap

        if self.kind == 'mass':
            constituents = self.config.property_package.component_list
            mass_basis = self.basis * (pyunits.kg/pyunits.hour)
            mass_in = 0
            for constituent in constituents:
                mass_in += self.conc_mass_in[time, constituent]
            density = 0.6312 * mass_in + 997.86
            total_mass_in = density * self.flow_in
            mass_factor = total_mass_in / mass_basis
            basic_cap = self.cap_basis * mass_factor ** self.cap_exp
            return basic_cap

    def elect(self):
        '''
        Electricity intensity for basic units.

        :return: Electricity intensity [kWh/m3]
        '''
        if self.unit_process_name in ['mbr_nitrification', 'mbr_denitrification'] and not self.case_specific:
            # Electricity consumption for MBRs from:
            # "Assessing Location and Scale of Urban Nonpotable Water Reuse Systems for Life-Cycle Energy Consumption and Greenhouse Gas Emissions" Kavvada et al (2016)
            # Equation located in SI
            return 9.5 * self.flow_in ** -0.3
        else:
            return self.elect_intensity

    def get_costing(self):
        '''
        Initialize the unit in WaterTAP3.
        '''
        time = self.flowsheet().config.time.first()
        self.flow_in = pyunits.convert(self.flow_vol_in[time],
            to_units=pyunits.m**3/pyunits.hr)
        self.unit_process_name = self.unit_params['unit_process_name']
        if 'case_specific' in self.unit_params.keys():
            self.case_specific = self.unit_params['case_specific']
            self.basis, self.cap_basis, self.cap_exp, \
            self.elect_intensity, self.basis_year, self.kind = cost_curves.basic_unit(
                self.unit_process_name, case_specific=self.case_specific)
        else:
            self.case_specific = False
            self.basis, self.cap_basis, self.cap_exp, \
            self.elect_intensity, self.basis_year, self.kind = cost_curves.basic_unit(self.unit_process_name)
        
        self.deltaP_outlet.unfix()
        self.deltaP_waste.unfix()
        self.pressure_out.fix(1)
        self.pressure_waste.fix(1)

        self.costing.fixed_cap_inv_unadjusted = Expression(expr=self.fixed_cap(),
                doc='Unadjusted fixed capital investment')
        self.electricity = Expression(expr=self.elect(),
                doc='Electricity intensity [kWh/m3]')
        financials.get_complete_costing(self.costing, basis_year=self.basis_year)

