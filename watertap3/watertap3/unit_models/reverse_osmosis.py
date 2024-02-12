from pyomo.environ import (
    Param,
    Block,
    Constraint,
    Expression,
    NonNegativeReals,
    check_optimal_termination,
    Var,
    units as pyunits,
)
from copy import deepcopy
from idaes.core.solvers.get_solver import get_solver
from idaes.core import declare_process_block_class
from watertap3.utils import financials
from watertap3.core.wt3_unit_sido import WT3UnitProcessSIDOData
from pyomo.network import Port
import idaes.logger as idaeslog
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

_log = idaeslog.getLogger(__name__)

module_name = "reverse_osmosis"

def cost_reverse_osmosis(blk):

    blk.membrane_cost = Var(
        initialize=40, bounds=(10, 80), domain=NonNegativeReals, doc="Membrane cost"
    )

    blk.factor_membrane_replacement = Var(
        initialize=0.2,
        domain=NonNegativeReals,
        bounds=(0.01, 3),
        doc="replacement rate membrane fraction",
    )
    blk.pressure_vessel_cost = Var(
        # initialize=2e6,
        bounds=(0, None),
        doc="Pressure vessel cost",
    )
    blk.rack_support_cost = Var(
        # initialize=1e6,
        bounds=(0, None),
        doc="Rack support cost",
    )

    blk.fix_all_vars()


@declare_process_block_class("UnitProcess")
class UnitProcessData(WT3UnitProcessSIDOData):
    def build(self):

        super().build()
        # self.del_component(self.properties_in)
        # self.feed = Block()
        # tmp_dict = dict(**self.config.property_package_args)
        # tmp_dict["has_phase_equilibrium"] = False
        # tmp_dict["parameters"] = self.config.property_package
        # tmp_dict["defined_state"] = True
        # self.feed.properties_in = prop_in = self.config.property_package.state_block_class(
        #     doc="Material properties of inlet stream", **tmp_dict
        # )
        # self.inlet.remove("pressure")
        # for arc in self.inlet.arcs():
        #     for port in arc.ports:
        #         port.remove("pressure")
        # self.del_component(self.inlet)
        # for port in self.component_objects(Port):
        #     for arc in port.arcs():
        #         for p in arc:
        #             print(port, arc, p)
        #             p.remove("pressure")
        #             p.remove("flow_vol")
        #             p.remove("conc_mass")


        self.pressure_atm = Param(
            initialize=101325,
            mutable=True,
            units=pyunits.Pa,
            doc="Atmospheric pressure",
        )
        self.module_membrane_area = Param(
            initialize=37.16,
            mutable=True,
            units=pyunits.m**3,
            doc="Membrane area per module",
        )

        self.water_permeability = Param(
            initialize=4.2e-12,
            mutable=True,
            units=pyunits.m / (pyunits.Pa * pyunits.s),
            doc="Water permeability",
        )

        self.salt_permeability = Param(
            initialize=3.5e-8,
            mutable=True,
            units=pyunits.m / pyunits.s,
            doc="Salt permeability",
        )

        # self.water_permeability = Var(
        #     initialize=4.2e-12,
        #     bounds=(1e-12, 9e-12),
        #     domain=NonNegativeReals,
        #     units=pyunits.m / (pyunits.Pa * pyunits.s),
        #     doc="Water permeability",
        # )

        # self.salt_permeability = Var(
        #     initialize=3.5e-8,
        #     bounds=(1e-8, 9e-8),
        #     domain=NonNegativeReals,
        #     units=pyunits.m / pyunits.s,
        #     doc="Salt permeability",
        # )
        self.deltaP_outlet = Var(
            initialize=1e-6,
            units=pyunits.Pa,
            doc="Pressure change between inlet and outlet",
        )

        # self.deltaP_outlet.fix(0)

        self.deltaP_waste = Param(
            initialize=3e5,
            units=pyunits.Pa,
            doc="Pressure change between inlet and waste",
        )
        # self.deltaP_waste.fix(0)
        self.membrane_area = Var(
            initialize=1e6,
            domain=NonNegativeReals,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Total membrane area",
        )

        # self.pressure_operating = Var(
        #     initialize=100,
        #     domain=NonNegativeReals,
        #     bounds=(0, None),
        #     units=pyunits.Pa,
        #     doc="Total membrane area",
        # )
        # self.deltaP_outlet.unfix()
        # self.deltaP_waste.fix(3e6)

        @self.Expression(doc="Number membrane modules needed")
        def number_modules(b):
            return b.membrane_area / b.module_membrane_area
        
        @self.Expression(doc="")
        def operating_pressure(b):
            return (
                (b.properties_in.pressure + b.properties_waste.pressure) * 0.5
            ) 

        @self.Constraint(doc="Permeate pressure")
        def eq_permeate_pressure(b):
            return b.deltaP_outlet == -b.properties_in.pressure + b.pressure_atm

        @self.Constraint(doc="Water transport equation")
        def eq_water_transport(b):
            avg_osm_pressure = (
                b.properties_in.pressure_osmotic + b.properties_waste.pressure_osmotic
            ) * 0.5
            # b.avg_pressure_in = 
            # avg_pressure_in = (
            #     (b.properties_in.pressure + b.properties_waste.pressure) * 0.5
            # ) 
            return (
                b.properties_out.flow_mass_comp["H2O"]
                == b.properties_in.dens_mass
                * b.water_permeability
                * (b.operating_pressure - avg_osm_pressure)
                * b.membrane_area
            )

        
        @self.Constraint(doc="Salt transport equation")
        def eq_salt_transport(b):
            return (
                b.properties_out.flow_mass_comp["tds"]
                == b.salt_permeability
                * (
                    b.properties_in.conc_mass_comp["tds"]
                    - b.properties_out.conc_mass_comp["tds"]
                )
                * b.membrane_area
            )

        @self.Constraint(doc="Pressure of brine stream")
        def eq_pressure_waste(b):
            return b.properties_waste.pressure == b.properties_in.pressure + b.deltaP_waste
        
        @self.Constraint(doc="Pressure of perm stream")
        def eq_pressure_permeate(b):
            return b.properties_out.pressure == b.properties_in.pressure + b.deltaP_outlet
        
        @self.Constraint()
        def eq_min_pressure(b):
            return b.properties_in.pressure >= b.properties_in.pressure_osmotic
        @self.Constraint()
        def eq_min_pressure_brine(b):
            return b.properties_waste.pressure >= b.properties_waste.pressure_osmotic
        
        # @self.Constraint(doc="Water recovery")
        # def eq_water_recovery(b):
        #     return (
        #         b.water_recovery == b.properties_out.flow_vol / b.properties_in.flow_vol
        #     )                     doc="Estimated density of solution")

    # def fixed_cap(self, t, b_cost):
    #     '''

    #     :param t: Indexing variable for Pyomo Var()
    #     :type t: int
    #     :param b_cost: Costing block for unit.
    #     :type b_cost: object
    #     :return: Fixed capital costs for reverse osmosis [$MM]
    #     '''
    #     ro_cap = (self.tpec_tic * (b_cost.pump_capital_cost + b_cost.mem_capital_cost + b_cost.erd_capital_cost) +
    #               3.3 * (self.pressure_vessel_cost[t] + self.rack_support_cost[t])) * 1E-6  # $MM ### 1.65 is TIC
    #     return ro_cap

    # def elect(self, t):
    #     '''

    #     :param t: Indexing variable for Pyomo Var()
    #     :type t: int
    #     :return:
    #     '''

    #     electricity = ((self.pump_power - self.erd_power) / 1000) / (self.flow_vol_in[t] * 3600)
    #     return electricity

    # def _set_pressure(self, b, ):
    #     b.pressure = Var(
    #                      initialize=45,
    #                      domain=NonNegativeReals,
    #                      bounds=(2, 90),
    #                      doc='pressure')

    # def _set_constraints(self, t):

    #     ## FLUX CONSTRAINTS
    #     self.water_salt_perm_eq1 = Constraint(
    #             expr=self.b[t] <= (0.083 * self.a[t] - 0.002) * 1.25)
    #     self.water_salt_perm_eq2 = Constraint(
    #             expr=self.b[t] >= (0.083 * self.a[t] - 0.002) * 0.75)

    #     ## RETENTATE CONSTRAINTS

    #     self.retentate.eq2 = Constraint(
    #             expr=self.retentate.conc_mass_H2O[t] == self.retentate.conc_mass_total[t] - self.conc_mass_waste[t, 'tds'])
    #     self.retentate.eq3 = Constraint(
    #             expr=self.retentate.mass_frac_tds[t] * self.retentate.conc_mass_total[t] == self.conc_mass_waste[t, 'tds'])

    #     ## PERMEATE CONSTRAINTS
    #     self.permeate.eq1 = Constraint(
    #             expr=self.permeate.conc_mass_total[t] == 756 * self.permeate.mass_frac_tds[t] * 1E-6 + 995)
    #     self.permeate.eq2 = Constraint(
    #             expr=self.conc_mass_out[t, 'tds'] == self.permeate.conc_mass_total[t] * self.permeate.mass_frac_tds[t] * 1E-6)
    #     self.permeate.eq3 = Constraint(
    #             expr=self.permeate.mass_flow_H2O[t] == self.membrane_area[t] * self.pure_water_flux[t])
    #     self.permeate.eq4 = Constraint(
    #             expr=self.permeate.mass_flow_tds[t] == 0.5 * self.membrane_area[t] * self.b[t] * 1E-7 * (self.conc_mass_in[t, 'tds'] + self.conc_mass_waste[t, 'tds']))
    #     self.permeate.eq5 = Constraint(
    #             expr=self.permeate.mass_frac_tds[t] * (self.permeate.mass_flow_tds[t] + self.permeate.mass_flow_H2O[t]) == 1E6 * self.permeate.mass_flow_tds[t])
    #     self.permeate.eq6 = Constraint(
    #             expr=self.pure_water_flux[t] == self.pw * self.a[t] * 1E-7 * ((self.feed.pressure[t] - self.p_atm - self.pressure_drop * 0.5) -
    #                     (self.feed.pressure_osm[t] + self.retentate.pressure_osm[t]) * 0.5))

    #     # PRESSURE BALANCE
    #     self.momentum_balance_eq = Constraint(
    #             expr=self.retentate.pressure[t] == self.feed.pressure[t] - self.pressure_drop)
    #     self.flow_vol_eq1 = Constraint(
    #             expr=self.flow_vol_out[t] * self.permeate.conc_mass_total[t] == (self.permeate.mass_flow_tds[t] + self.permeate.mass_flow_H2O[t]))
    #     self.flow_vol_eq2 = Constraint(
    #             expr=self.flow_vol_waste[t] * self.retentate.conc_mass_total[t] == (self.retentate.mass_flow_tds[t] + self.retentate.mass_flow_H2O[t]))
    #     self.pressure_waste_outlet_eq = Constraint(
    #             expr=self.feed.pressure[t] - self.pressure_drop == self.pressure_waste[t])

    #     # MASS BALANCE
    #     self.mass_balance_H2O = Constraint(
    #             expr=self.feed.mass_flow_H2O[t] == self.permeate.mass_flow_H2O[t] + self.retentate.mass_flow_H2O[t])
    #     self.mass_balance_tds = Constraint(
    #             expr=self.feed.mass_flow_tds[t] == self.permeate.mass_flow_tds[t] + self.retentate.mass_flow_tds[t])

    #     # PERMEATE PRESSURE
    #     self.p_out_eq = Constraint(
    #             expr=1 == self.pressure_out[t])
    #     self.pump_constraint_power = Constraint(
    #             expr=self.pump_power >= 0)

    #     # VESSEL COST
    #     self.pressure_vessel_cost_eq = Constraint(
    #             expr=self.pressure_vessel_cost[t] * 0.99 <= self.membrane_area[t] * 0.025 * 1000)  # assumes 2 trains. 150 ft start, 5ft per additional vessel. EPA.
    #     self.rack_support_cost_eq = Constraint(
    #             expr=self.rack_support_cost[t] * 0.99 <= (150 + (self.membrane_area[t] * 0.025 * 5)) * 33 * 2)
    #     self.pressure_vessel_cost_eq2 = Constraint(
    #             expr=self.pressure_vessel_cost[t] * 1.01 >= self.membrane_area[t] * 0.025 * 1000)  # assumes 2 trains. 150 ft start, 5ft per additional vessel. EPA.
    #     self.rack_support_cost_eq2 = Constraint(
    #             expr=self.rack_support_cost[t] * 1.01 >= (150 + (self.membrane_area[t] * 0.025 * 5)) * 33 * 2)

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for initialization routines

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        flags = self.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        # init_log.info("Initialization Step 1a Complete.")

        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # Set state_args from inlet state
        if state_args is None:
            self.state_args = state_args = {}
            state_dict = self.properties_in.define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value


        self.state_args_out = state_args_out = deepcopy(state_args)
        for k, v  in state_args.items():
            if k == "flow_vol":
                state_args_out[k] == v * 0.5
            elif k == "conc_mass_comp":
            # elif isinstance(v, dict):
                state_args_out[k] == dict()
                for j, u in v.items():
                    state_args_out[k][j] = (1 - self.removal_fraction[j].value) * u
                    # if j == "tds":
                    #     self.removal_fraction[j].unfix()
                        
                        

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        self.state_args_waste = state_args_waste = deepcopy(state_args)
        for k, v  in state_args.items():
            if k == "flow_vol":
                state_args_waste[k] == v * 0.5
            elif k == "conc_mass_comp":
            # elif isinstance(v, dict):
                state_args_waste[k] == dict()
                for j, u in v.items():
                    state_args_waste[k][j] = self.removal_fraction[j].value * u
        
        self.properties_waste.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_waste,
        )

        init_log.info("Initialization Step 1c Complete.")

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)

        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Release Inlet state
        self.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        self.removal_fraction["tds"].unfix()

        # if not check_optimal_termination(res):
        #     raise InitializationError(f"Unit model {self.name} failed to initialize.")


    def calculate_scaling_factors(self):
        print("\n\n\ncalc scaling factors RO\n\n")

    @property
    def default_costing_method(self):
        return cost_reverse_osmosis

    def get_costing(self):
        """
        Initialize the unit in WaterTAP3.
        """

        basis_year = 2007
        tpec_tic = "TIC"

        sys_cost_params = self.parent_block().costing_param
        self.parent_block().has_ro = True
        if "erd" in self.unit_params.keys():
            self.erd = self.unit_params["erd"]
            if self.erd not in ["yes", "no"]:
                self.erd = "no"
        else:
            self.erd = "no"
        # self.del_component(self.outlet_pressure_constraint)
        # self.del_component(self.waste_pressure_constraint)
        # self.del_component(self.recovery_equation)
        # self.del_component(self.flow_balance)
        # self.del_component(self.component_removal_equation)

        self.deltaP_waste.unfix()
        self.deltaP_outlet.unfix()

        self.permeate = Block()
        self.feed = Block()
        self.retentate = Block()

        self.units_meta = self.config.property_package.get_metadata().get_derived_units

        # DEFINE VARIABLES
        self.feed.water_flux = Var(
            initialize=5e-3,
            bounds=(1e-5, 1.5e-2),
            units=self.units_meta("mass")
            * self.units_meta("length") ** -2
            * self.units_meta("") ** -1,
            domain=NonNegativeReals,
            doc="water flux",
        )
        self.retentate.water_flux = Var(
            initialize=5e-3,
            bounds=(1e-5, 1.5e-2),
            units=self.units_meta("mass")
            * self.units_meta("length") ** -2
            * self.units_meta("") ** -1,
            domain=NonNegativeReals,
            doc="water flux",
        )
        self.pure_water_flux = Var(
            initialize=5e-3,
            bounds=(1e-3, 1.5e-2),
            units=self.units_meta("mass")
            * self.units_meta("length") ** -2
            * self.units_meta("") ** -1,
            domain=NonNegativeReals,
            doc="water flux",
        )
        self.a = Var(
            initialize=4.2,
            bounds=(1, 9),
            domain=NonNegativeReals,
            doc="water permeability",
        )
        self.b = Var(
            initialize=0.35,
            bounds=(0.1, 0.9),
            domain=NonNegativeReals,
            doc="Salt permeability",
        )
        self.mem_cost = Var(
            initialize=40, bounds=(10, 80), domain=NonNegativeReals, doc="Membrane cost"
        )
        self.membrane_area = Var(
            initialize=1e5, domain=NonNegativeReals, bounds=(1e1, 1e12), doc="area"
        )
        self.factor_membrane_replacement = Var(
            initialize=0.2,
            domain=NonNegativeReals,
            bounds=(0.01, 3),
            doc="replacement rate membrane fraction",
        )
        self.pressure_vessel_cost = Var(
            # initialize=2e6,
            bounds=(0, None),
            doc="Pressure vessel cost",
        )
        self.rack_support_cost = Var(
            # initialize=1e6,
            bounds=(0, None),
            doc="Rack support cost",
        )
        self.factor_membrane_replacement.fix(0.25)  #

        # from excel regression based on paper for membrane cost y = 0.1285x - 0.0452 #R² =
        # 0.9932. y = b. x = a.
        # same 4 membranes used for regression ŷ = 15.04895X1 - 131.08641X2 + 29.43797

        for b in [self.permeate, self.feed, self.retentate]:
            self._set_flow_mass(
                b,
            )
            self._set_mass_frac(
                b,
            )
            self._set_conc_mass(
                b,
            )
            self._set_osm_coeff(
                b,
            )
            self._set_pressure_osm(
                b,
            )
            if str(b) == "permeate":
                continue
            else:
                self._set_pressure(
                    b,
                )

        self.feed.pressure.unfix()
        self.mem_cost.fix(30)


        ## CONSTANTS
        self.pump_eff = 0.8
        self.erd_eff = 0.9
        self.pressure_drop = 3 * pyunits.bar
        self.p_atm = 1
        self.pw = 1000

        self.pressure_diff = (
            self.feed.pressure[t] - self.pressure_in[t]
        ) * 1e5  # assumes atm pressure before pump. change to Pa
        self.pump_power = (self.flow_vol_in[t] * self.pressure_diff) / self.pump_eff

        self._set_constraints(t)

        flow_in_m3hr = pyunits.convert(
            self.flow_vol_in[t], to_units=pyunits.m**3 / pyunits.hour
        )
        flow_out_m3hr = pyunits.convert(
            self.flow_vol_out[t], to_units=pyunits.m**3 / pyunits.hour
        )
        flow_waste_m3hr = pyunits.convert(
            self.flow_vol_waste[t], to_units=pyunits.m**3 / pyunits.hour
        )

        try:
            scaling = self.unit_params["scaling"]

        except:
            scaling = "no"

        for j in self.const_list2:
            if scaling == "yes":
                self.del_component(self.component_mass_balance)
                # setattr(self, ('%s_eq1' % j), Constraint(
                #         expr=self.flow_vol_in[t] * self.conc_mass_in[j] == self.flow_vol_out[t] *
                #              self.conc_mass_out[j] + self.flow_vol_waste[t] * self.conc_mass_waste[j]))

                # setattr(self, ('%s_eq' % j), Constraint(
                #         expr=self.removal_fraction[j] * self.flow_vol_in[t] * self.conc_mass_in[j] == self.flow_vol_waste[t] * self.conc_mass_waste[j]))
                setattr(
                    self,
                    ("%s_eq1" % j),
                    Constraint(
                        expr=flow_in_m3hr
                        * pyunits.convert(
                            self.conc_mass_in[j], to_units=pyunits.mg / pyunits.liter
                        )
                        == flow_out_m3hr
                        * pyunits.convert(
                            self.conc_mass_out[j], to_units=pyunits.mg / pyunits.liter
                        )
                        + flow_waste_m3hr
                        * pyunits.convert(
                            self.conc_mass_out[j], to_units=pyunits.mg / pyunits.liter
                        )
                    ),
                )
                setattr(
                    self,
                    ("%s_eq" % j),
                    Constraint(
                        expr=self.removal_fraction[j]
                        * flow_in_m3hr
                        * pyunits.convert(
                            self.conc_mass_in[j], to_units=pyunits.mg / pyunits.liter
                        )
                        == flow_waste_m3hr
                        * pyunits.convert(
                            self.conc_mass_waste[j], to_units=pyunits.mg / pyunits.liter
                        )
                    ),
                )
            else:
                setattr(
                    self,
                    ("%s_eq" % j),
                    Constraint(
                        expr=self.removal_fraction[j]
                        * flow_in_m3hr
                        * pyunits.convert(
                            self.conc_mass_in[j], to_units=pyunits.mg / pyunits.liter
                        )
                        == flow_waste_m3hr
                        * pyunits.convert(
                            self.conc_mass_waste[j], to_units=pyunits.mg / pyunits.liter
                        )
                    ),
                )

        b_cost = self.costing
        b_cost.pump_capital_cost = self.pump_power * (53 / 1e5 * 3600) ** 0.97
        b_cost.pressure_vessel_cap_cost = (
            self.pressure_vessel_cost[t] + self.rack_support_cost[t]
        )

        ################ Energy Recovery
        # assumes atmospheric pressure out
        if self.erd == "yes":
            x_value = (
                (self.retentate.mass_flow_tds[t] + self.retentate.mass_flow_H2O[t])
                / self.retentate.conc_mass_total[t]
                * 3600
            )
            b_cost.erd_capital_cost = (
                3134.7 * (x_value * self.retentate.conc_mass_total[t]) ** 0.58
            )
            self.erd_power = (
                self.flow_vol_waste[t] * (self.retentate.pressure[t] - 1) * 1e5
            ) / self.erd_eff

        if self.erd == "no":
            self.erd_power = 0

        b_cost.erd_capital_cost = 0
        b_cost.mem_capital_cost = self.mem_cost[t] * self.membrane_area[t]

        self.costing.fixed_cap_inv_unadjusted = Expression(
            expr=self.fixed_cap(t, b_cost), doc="Unadjusted fixed capital investment"
        )

        self.num_membranes = (
            self.membrane_area[t] / 37.161
        )  # 400 ft2 / membrane = 37.161 m2 / membrane from BRACKISH paper
        self.ro_recovery = self.flow_vol_out[t] / self.flow_vol_in[t]
        self.flux = self.flow_vol_out[t] / (self.membrane_area[t] * pyunits.m**2)
        self.flux_lmh = pyunits.convert(
            self.flux, to_units=(pyunits.liter / pyunits.m**2 / pyunits.hour)
        )
        self.salt_rejection_mass = (
            1 - self.permeate.mass_flow_tds[t] / self.feed.mass_flow_tds[t]
        ) * 100
        self.salt_rejection_conc = (
            1 - self.conc_mass_out[t, "tds"] / self.conc_mass_in[t, "tds"]
        ) * 100
        self.SEC = (
            (self.flow_vol_in[t] * self.pressure_drop) / self.flow_vol_out[t]
        ) * self.pump_eff  # specific energy consumption, from - Zhu, A., Christofides, P., Cohen, Y., "On RO membrane and energy costs..."
        self.salt_flux = self.b[t] * (
            self.conc_mass_in[t, "tds"] - self.conc_mass_out[t, "tds"]
        )

        ################ operating
        # membrane operating cost
        b_cost.other_var_cost = (
            self.factor_membrane_replacement[t]
            * self.mem_cost[t]
            * self.membrane_area[t]
            * sys_cost_params.plant_cap_utilization
            * 1e-6
        )
        self.electricity = Expression(
            expr=self.elect(t), doc="Electricity intensity [kWh/m3]"
        )
        ####### electricity and chems
        sys_specs = self.parent_block().costing_param
        # self.electricity = ((self.pump_power - self.erd_power) / 1000) / (self.flow_vol_in[t] * 3600)
        b_cost.pump_electricity_cost = (
            1e-6 * (self.pump_power / 1000) * 365 * 24 * sys_specs.electricity_price
        )
        b_cost.erd_electricity_sold = (
            1e-6 * (self.erd_power / 1000) * 365 * 24 * sys_specs.electricity_price
        )
        # b_cost.electricity_cost = (b_cost.pump_electricity_cost - b_cost.erd_electricity_sold) * sys_cost_params.plant_cap_utilization

        self.chem_dict = {"unit_cost": 0.01}

        financials.get_complete_costing(
            self.costing, basis_year=basis_year, tpec_tic=tpec_tic
        )
