import ast
import os
from copy import deepcopy

import numpy as np
import pandas as pd

from pyomo.environ import Block, ConcreteModel, Var
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import *
from . import module_import

# from .source_wt3 import Source
from ..core.source_wt3 import Source
from ..core.mixer_wt3 import Mixer
from ..core.splitter_wt3 import Splitter
from ..core.wt3_property_package import WT3ParameterBlock
from watertap3.core.wt3_unit_pt import WT3UnitProcessPT
from watertap3.core.wt3_unit_siso import WT3UnitProcessSISOData
from watertap3.core.wt3_unit_sido import WT3UnitProcessSIDOData
from watertap3.unit_models.reverse_osmosis import ReverseOsmosis
from ..costing import WT3Costing, WT3UnitCosting

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

treatment_train_setup_file = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/data/treatment_train_setup.csv"
)
case_study_water_sources_file = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/data/case_study_water_sources.csv"
)
water_recovery_file = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/data/water_recovery_factors.csv"
)
wr_df = pd.read_csv(water_recovery_file)
constituent_removal_factors_file = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/data/constituent_removal_factors.csv"
)
crf_df = pd.read_csv(constituent_removal_factors_file)

__all__ = [
    "watertap3_setup",
    "get_case_study",
    "add_unit_process",
    "add_water_source",
    "get_pfd_dict",
    #    'generate_constituent_list',
    "create_wt3_unit",
    "create_arcs",
    "create_arc_dict",
    "check_split_mixer_need",
    "create_mixers",
    "create_splitters",
    "add_waste_streams",
]


def watertap3_setup(
    case_study=None,
    reference="nawi",
    scenario="baseline",
    source_reference=None,
    source_case_study=None,
    source_scenario=None,
    print_it=True,
    ro_bounds="swro",
    desired_recovery=1,
):
    """
    Initial setup of WaterTAP3 model.

    Create flowsheet and read in basic information about model (water sources, units in treatment train)
    """

    def get_source(reference, water_type, case_study, scenario):
        """
        Read in water source data.
        """
        df = pd.read_csv(case_study_water_sources_file, index_col="variable")
        if "test" in case_study:
            case_study = "test"
        try:
            source_df = df[
                (
                    (df.case_study == case_study)
                    & (df.water_type == water_type)
                    & (df.reference == reference)
                    & (df.scenario == scenario)
                )
            ].copy()
            source_flow = source_df.loc["flow"].value
        except:
            source_df = df[
                (
                    (df.case_study == case_study)
                    & (df.water_type == water_type)
                    & (df.reference == reference)
                    & (df.scenario == "baseline")
                )
            ].copy()
            source_flow = source_df.loc["flow"].value
        # source_df.drop(source_df[source_df.index == "flow"].index, inplace=True)
        return source_flow, source_df

    case_study_print = case_study.replace("_", " ").swapcase()
    scenario_print = scenario.replace("_", " ").swapcase()

    m_name = f"{case_study_print}: {scenario_print}"
    m = ConcreteModel(name=m_name)

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.train = {
        "case_study": case_study,
        "reference": reference,
        "scenario": scenario,
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

    df = pd.read_csv(treatment_train_setup_file)  # Read in treatment train input sheet.

    water_type_list = []
    # if new_df_units is not None:
    #     m.fs.df_units = new_df_units.copy()
    # else:
    m.fs.df_units = df[
        (
            (df.Reference == reference)
            & (df.Scenario == scenario)
            & (df.CaseStudy == case_study)
        )
    ].copy()

    # if 'binary_train' in m.fs.df_units.columns:
    #     cols = ['CaseStudy', 'binary_train', 'Reference', 'Scenario', 'Unit', \
    #         'Type', 'UnitName', 'ToUnitName', 'FromPort', 'Parameter']
    # else:
    cols = [
        "CaseStudy",
        "Reference",
        "Scenario",
        "Unit",
        "Type",
        "UnitName",
        "ToUnitName",
        "FromPort",
        "Parameter",
    ]

    m.fs.df_units = m.fs.df_units[cols].copy()

    m.fs.has_ro = False
    m.fs.has_ix = False
    if "ion_exchange" in m.fs.df_units.Unit:
        m.fs.has_ix = True
    if "reverse_osmosis" in m.fs.df_units.Unit:
        m.fs.has_ro = True

    for i in m.fs.df_units[m.fs.df_units.Type == "intake"].index:
        temp_dict = ast.literal_eval(
            m.fs.df_units[m.fs.df_units.Type == "intake"].loc[i]["Parameter"]
        )
        for water_type in temp_dict["water_type"]:
            water_type_list.append(water_type)

    if len(water_type_list) == 1:
        water_type_list = water_type_list[0]

    m.fs.source_water = {
        "case_study": source_case_study,
        "reference": source_reference,
        "scenario": source_scenario,
        "water_type": water_type_list,
    }

    flow_dict = {}

    if isinstance(m.fs.source_water["water_type"], list):
        # m.fs.source_df = pd.DataFrame()
        tmp = list()
        for water_type in m.fs.source_water["water_type"]:
            source_flow, source_df = get_source(
                m.fs.source_water["reference"],
                water_type,
                m.fs.source_water["case_study"],
                m.fs.source_water["scenario"],
            )
            flow_dict[water_type] = source_flow
            # m.fs.source_df = m.fs.source_df.append(source_df)
            tmp.append(source_df)
        m.fs.source_df = pd.concat(tmp)
    else:
        source_flow, source_df = get_source(
            m.fs.source_water["reference"],
            m.fs.source_water["water_type"],
            m.fs.source_water["case_study"],
            m.fs.source_water["scenario"],
        )
        flow_dict[m.fs.source_water["water_type"]] = source_flow
        m.fs.source_df = source_df

    if print_it:
        print(f"\nCase Study = {case_study_print}" f"\nScenario = {scenario_print}")

    m.fs.flow_in_dict = flow_dict
    m.fs.unit_removal_fractions = dict()
    m.fs.unit_water_recovery = dict()
    x = generate_constituent_list(m)
    # x += ["H2O"]
    m.fs.water = WT3ParameterBlock(constituent_list=x)
    m.fs.costing = WT3Costing()

    return m


def get_case_study(m, new_df_units=None, print_it=True):
    """
    Function to add constituents and unit processes to flowsheet and connect
    all ports.
    """

    m.fs.pfd_dict = pfd_dict = get_pfd_dict(m.fs.df_units)

    # if print_it:
    print("\n=========================ADDING UNIT PROCESSES=========================")
    for unit_process_name in pfd_dict.keys():
        unit = unit_process_name.replace("_", " ").swapcase()
        unit_process_type = pfd_dict[unit_process_name]["Unit"]
        unit_process_kind = pfd_dict[unit_process_name]["Type"]
        print(f"{unit}")
        m = add_unit_process(
            m=m,
            unit_process_name=unit_process_name,
            unit_process_type=unit_process_type,
            unit_process_kind=unit_process_kind,
        )
        # u = getattr(m.fs, unit_process_name)
        # u.properties_in.flow_vol
        # u.properties_in.conc_mass_comp[...]
    print("=======================================================================\n")
    # else:
    #     for unit_process_name in pfd_dict.keys():
    #         unit = unit_process_name.replace("_", " ").swapcase()
    #         unit_process_type = pfd_dict[unit_process_name]["Unit"]
    #         unit_process_kind = pfd_dict[unit_process_name]["Type"]
    #         m = add_unit_process(
    #             m=m,
    #             unit_process_name=unit_process_name,
    #             unit_process_type=unit_process_type,
    #             unit_process_kind=unit_process_kind,
    #         )

    # create a dictionary with all the arcs in the network based on the pfd_dict
    m, arc_dict, arc_i = create_arc_dict(m, pfd_dict, m.fs.flow_in_dict)
    m.fs.arc_dict = arc_dict
    # calculate_scaling_factors(m)
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


def add_unit_process(
    m=None, unit_process_name=None, unit_process_type=None, unit_process_kind=None
):
    # print(unit_process_type)
    UnitProcess = module_import.get_module(unit_process_type)

    unit_params = m.fs.pfd_dict[unit_process_name]["Parameter"]
    if not isinstance(unit_params, dict):
        unit_params = dict()

    if unit_process_type == "basic_unit":
        setattr(
            m.fs,
            unit_process_name,
            UnitProcess(property_package=m.fs.water, unit_params=unit_params),
        )
        basic_unit_name = unit_params["unit_process_name"]
        m = create_wt3_unit(m, basic_unit_name, unit_process_name)

    else:
        setattr(
            m.fs,
            unit_process_name,
            UnitProcess(property_package=m.fs.water, unit_params=unit_params),
        )
        # if not isinstance(getattr(m.fs, unit_process_name), WT3UnitProcessPT):
        m = create_wt3_unit(m, unit_process_type, unit_process_name)

    unit = getattr(m.fs, unit_process_name)
    # unit.chem_dict = {}

    # unit.unit_params = unit_params
    unit.costing = WT3UnitCosting(flowsheet_costing_block=m.fs.costing)
    # unit.costing.get_cost_indices()

    unit.unit_type = unit_process_type
    unit.unit_name = unit_process_name
    unit.unit_kind = unit_process_kind
    # TODO: set this in unit model file rather than here
    unit.unit_pretty_name = (
        unit_process_name.replace("_", " ")
        .title()
        .replace("Ro", "RO")
        .replace("Zld", "ZLD")
        .replace("Aop", "AOP")
        .replace("Uv", "UV")
        .replace("And", "&")
        .replace("Sw", "SW")
        .replace("Gac", "GAC")
        .replace("Ph", "pH")
        .replace("Bc", "BC")
        .replace("Wwtp", "WWTP")
        .replace("Pac", "PAC")
        .replace("Co2", "CO2")
        .replace("Kmno4", "KMnO4")
    )

    return m


def add_water_source(m=None, source_name=None, water_type=None, flow=None):
    setattr(m.fs, source_name, Source(property_package=m.fs.water))

    source = getattr(m.fs, source_name)
    # source.properties.flow_vol.fix(flow)
    
    # source.properties.flow_mass_comp[...]
    temp_source_df = m.fs.source_df[m.fs.source_df.water_type == water_type].copy()
    train_constituent_list = list()
    var_args = dict()
    for constituent_name in source.config.property_package.component_list:
        # print(constituent_name)
        # var_args[("conc_mass_comp", constituent_name)] = 0
        if constituent_name == "H2O":
            var_args[("flow_vol", None)] = flow
            source.properties.flow_mass_comp[constituent_name].set_value(1000 * flow)
            # source.properties.conc_mass_comp[constituent_name].set_value(1000)
        elif constituent_name in temp_source_df.index:
            conc = temp_source_df.loc[constituent_name].value
            var_args[("conc_mass_comp", constituent_name)] = conc
            # source.properties.conc_mass_comp[constituent_name].fix(conc)
            source.properties.flow_mass_comp[constituent_name].set_value(conc * flow)
        else:
            pass
            # source.properties.flow_mass_comp[constituent_name].fix(0)
            # source.properties.conc_mass_comp[constituent_name].fix(0)
    source.properties.flow_vol
    source.properties.conc_mass_comp[...]
    var_args[("pressure", None)] = 101325
    source.properties.calculate_state(var_args=var_args, hold_state=True)
    return m


def create_wt3_unit(m, unit_process_type, unit_process_name):
    def get_removal_factors():
        train = m.fs.train
        conditions = (
                (crf_df.unit_process == unit_process_type)
                & (crf_df.scenario == "baseline")
                & (crf_df.reference == train["reference"])
                & (crf_df.case_study.isin(["default", train["case_study"]]))
            )
        df = crf_df[conditions].copy()
        cl = getattr(m.fs, unit_process_name).config.property_package.solute_set
        removal_dict = dict()
        for c in cl:
            # removal factor df specific to this case study and constituent
            tmp = df[
                ((df.case_study == train["case_study"]) & (df.constituent == c))
            ].copy()
            if c not in df.constituent.unique(): 
                # if constituent isn't there, assume no removal
                rf = 0
            elif tmp.empty:
                # if there are no removal factors for this constituent specific to this case study,
                # use the default removal factor
                rf = df[
                    ((df.case_study == "default") & (df.constituent == c))
                ].value.iloc[0]
            else:
                # use case study specific removal factor for this unit
                rf = tmp.value.iloc[0]
            removal_dict[c] = float(rf)

        return removal_dict

    def get_water_recovery():

        conditions = (
            (wr_df.unit_process == unit_process_type)
            & (wr_df.case_study == m.fs.train["case_study"])
            & (wr_df.scenario == m.fs.train["scenario"])
        )
        unit_wr_df = wr_df[conditions].copy()
        if unit_wr_df.empty:
            # If not, we look for the default water recovery for this case study
            conditions = (
                (wr_df.unit_process == unit_process_type)
                & (wr_df.case_study == m.fs.train["case_study"])
                & (wr_df.scenario == "default")
            )
            unit_wr_df = wr_df[conditions].copy()
            if unit_wr_df.empty:
                # If not, we use the default water recovery for this unit
                conditions = (
                    (wr_df.unit_process == unit_process_type)
                    & (wr_df.case_study == "default")
                    & (wr_df.scenario == "default")
                )
                unit_wr_df = wr_df[conditions].copy()
                if unit_wr_df.empty:
                    raise ValueError(
                        f"No default water recovery factor provided for {unit_process_type}."
                        f"Please provide one in data/water_recovery_factors.csv."
                    )

        if len(unit_wr_df) == 1:
            pass

        elif len(unit_wr_df) > 1:
            if not len(unit_wr_df.recovery.unique()) == 1:
                raise ValueError(
                    f"Can only specify one water recovery factor per unit in water_recovery.csv,\n"
                    f"\tbut found {len(unit_wr_df)} for {unit_process_type}."
                )
            else:
                # If there are two identical entries
                unit_wr_df = unit_wr_df.iloc[0]

        unit_wrf = unit_wr_df.recovery.iloc[0]
        try:
            unit_wrf = float(unit_wrf)
            return unit_wrf
        except ValueError:
            print(f"water recovery unfixed for {unit_process_name}")
            return None

    unit = getattr(m.fs, unit_process_name)

    if hasattr(unit, "water_recovery"):
        unit_wrf = get_water_recovery()
        if unit_wrf is not None:
            unit.water_recovery.fix(unit_wrf)
        m.fs.unit_water_recovery[unit_process_name] = unit_wrf

    if hasattr(unit, "removal_fraction"):
        removal_dict = get_removal_factors()
        m.fs.unit_removal_fractions[unit_process_name] = removal_dict
        for c, rf in removal_dict.items():
            unit.removal_fraction[c].fix(rf)

    return m


def generate_constituent_list(m):
    train = m.fs.train
    # getting the list of consituents with removal factors that are bigger than 0
    df = pd.read_csv(constituent_removal_factors_file)
    df.case_study = np.where(
        df.case_study == "default", train["case_study"], df.case_study
    )
    df = df[df.reference == train["reference"]]
    df = df[df.case_study == train["case_study"]]
    df = df[df.scenario == "baseline"]
    list1 = df[df.value >= 0].constituent.unique()
    list2 = m.fs.source_df.index.unique().to_list()

    m.fs.constituent_list = constituent_list = [x for x in list1 if x in list2]

    return constituent_list


def get_pfd_dict(df_units):
    ### create pfd_dictionary for treatment train
    pfd_dict = df_units.set_index("UnitName").T.to_dict()
    for key in pfd_dict.keys():
        # parameter from string to dict
        if pfd_dict[key]["Parameter"] is not np.nan:
            pfd_dict[key]["Parameter"] = ast.literal_eval(pfd_dict[key]["Parameter"])
        if pfd_dict[key]["ToUnitName"] is not np.nan:
            if "," in pfd_dict[key]["ToUnitName"]:
                pfd_dict[key]["ToUnitName"] = pfd_dict[key]["ToUnitName"].split(",")
                pfd_dict[key]["FromPort"] = pfd_dict[key]["FromPort"].split(",")

    return pfd_dict


def create_arcs(m, arc_dict):
    for key in arc_dict.keys():
        source = arc_dict[key][0]
        source_port = arc_dict[key][1]
        outlet = arc_dict[key][2]
        outlet_port = arc_dict[key][3]
        setattr(
            m.fs,
            (f"arc{key}"),
            Arc(
                source=getattr(getattr(m.fs, source), source_port),
                destination=getattr(getattr(m.fs, outlet), outlet_port),
            ),
        )
        tmp = getattr(m.fs, f"arc{key}")
        tmp.propagated = False
    return m


def create_arc_dict(m, pfd_dict, flow):
    # create arc dictionary, add sources, add source to inlet arcs
    arc_i = 1
    arc_dict = dict()
    for u, d in pfd_dict.items():
        if d["Type"] == "intake":
            for water_type in d["Parameter"]["water_type"]:
                m = add_water_source(m=m, source_name=water_type, water_type=water_type, flow=flow[water_type])
                arc_dict[arc_i] = [water_type, "outlet", u, "inlet"]
                arc_i += 1
    for u, d in pfd_dict.items():
        if d["FromPort"] is not np.nan:
            if isinstance(d["FromPort"], list):
                for port_i in range(0, len(d["FromPort"])):
                    arc_dict[arc_i] = [u, d["FromPort"][port_i], d["ToUnitName"][port_i], "inlet"]
                    arc_i += 1
            else:
                arc_dict[arc_i] = [u, d["FromPort"], d["ToUnitName"], "inlet"]
                arc_i += 1

    return m, arc_dict, arc_i


def check_split_mixer_need(arc_dict):
    # check if a mixer or splitter is needed

    mixer_list = list() # (unit_name, inlet_port) pairs that need a mixer
    splitter_list = list() # (unit_name, outlet_port) pairs that need a splitter
    u1 = list() # unique (unit_name, inlet_port) pairs
    u2 = list() # unique (unit_name, outlet_port) pairs

    for (source, source_outlet, destination, destination_port) in arc_dict.values():
        if [source, source_outlet] not in u1:
            u1.append([source, source_outlet])
        elif [source, source_outlet] not in splitter_list:
                splitter_list.append([source, source_outlet])
        
        if [destination, destination_port] not in u2:
            u2.append([destination, destination_port])
        elif [destination, destination_port] not in mixer_list:
                mixer_list.append([destination, destination_port])
    
    return splitter_list, mixer_list


def create_mixers(m, mixer_list, arc_dict, arc_i):
    inlet_i = 1
    mixer_i = 1
    m.fs.mixers = list()
    for mixer_i, (unit, port) in enumerate(mixer_list, 1):
        inlet_list = list()
        mixer_name = f"mixer{mixer_i}"
        m.fs.mixers.append(mixer_name)
        # for k, v in arc_dict.items():
        for k in list(arc_dict.keys()):
            v = arc_dict[k]
            if (v[2] == unit) & (v[3] == port):
                inlet_name = f"inlet{inlet_i}"
                inlet_list.append(inlet_name)
                inlet_i += 1
                arc_dict[arc_i] = [v[0], v[1], mixer_name, inlet_name]
                arc_i += 1
                del arc_dict[k]

        # add mixer to model with inlet list
        mixer_config = dict(property_package=m.fs.water, inlet_list=inlet_list)
        setattr(m.fs, mixer_name, Mixer(**mixer_config))

        # arc from mixer outlet to node
        arc_dict[arc_i] = [mixer_name, "outlet", unit, port]
        arc_i += 1
        mixer_i += 1

    return m, arc_dict, mixer_i, arc_i


def create_splitters(m, splitter_list, arc_dict, arc_i):

    splitter_i = 1
    for splitter_i, (unit, port) in enumerate(splitter_list, 1):
        outlet_i = 1
        outlet_list = list()
        outlet_dict = dict()
        splitter_name = f"splitter{splitter_i}"
        for k in list(arc_dict.keys()):
            d = arc_dict[k]
            if (d[0] == unit) & (d[1] == port):
                split_dict = dict()
                w = 0
                for u in m.fs.pfd_dict[unit]["ToUnitName"]: 
                    if m.fs.pfd_dict[unit]["FromPort"][w] == "outlet":
                        if "split_fraction" in m.fs.pfd_dict[unit]["Parameter"]:
                            split_dict[u] = m.fs.pfd_dict[unit]["Parameter"]["split_fraction"][w]
                            w += 1
                outlet_name = f"outlet_{outlet_i}"
                outlet_list.append(outlet_name)
                outlet_i += 1
                arc_dict[arc_i] = [splitter_name, outlet_name, d[2], d[3]]
                arc_i += 1
                
                tmp = d[2]
                if d[2] not in split_dict.keys():
                    for x in arc_dict.keys():
                        y = arc_dict[x]
                        if d[2] == y[0]:
                            tmp = y[2]
                if len(split_dict) > 0:
                    outlet_dict[outlet_name] = split_dict[tmp]
                else:
                    outlet_dict[outlet_name] = "NA"
                del arc_dict[k]
        setattr(
            m.fs,
            splitter_name,
            Splitter(property_package=m.fs.water, outlet_dict=outlet_dict),
        )
        setattr(m.fs, f"{splitter_name}_outlet_dict", outlet_dict)
        arc_dict[arc_i] = [unit, port, splitter_name, "inlet"]
        arc_i += 1
        splitter_i += 1

    return m, arc_dict, splitter_i, arc_i


def add_waste_streams(m, arc_i, pfd_dict, mixer_i):
    # get number of units going to automatic waste disposal units
    i = 0
    m.fs.waste_inlet_list = waste_inlet_list = []
    unit_list = []
    for key in m.fs.pfd_dict.keys():
        unit_list.append(m.fs.pfd_dict[key]["Unit"])
        if "surface_discharge" == m.fs.pfd_dict[key]["Unit"]:
            sd_name = key
    if "surface_discharge" in unit_list:
        for b_unit in m.fs.component_objects(Block, descend_into=False):
            if hasattr(b_unit, "waste"):
                if len(getattr(b_unit, "waste").arcs()) == 0:
                    if str(b_unit)[3:] in list(pfd_dict.keys()):  #
                        if pfd_dict[str(b_unit)[3:]]["Type"] == "treatment":
                            i += 1
                            waste_inlet_list.append((f"inlet{i}"))

        if len(waste_inlet_list) > 1:
            i = 0
            waste_mixer = f"mixer{mixer_i}"
            m.fs.water_mixer_name = (
                waste_mixer  # used for displaying train. not used for model
            )
            setattr(
                m.fs,
                waste_mixer,
                Mixer(property_package=m.fs.water, inlet_list=waste_inlet_list),
            )
            for b_unit in m.fs.component_objects(Block, descend_into=False):
                if hasattr(b_unit, "waste"):
                    if len(getattr(b_unit, "waste").arcs()) == 0:
                        if str(b_unit)[3:] in list(pfd_dict.keys()):
                            if pfd_dict[str(b_unit)[3:]]["Type"] == "treatment":
                                i += 1
                                setattr(
                                    m.fs,
                                    (f"arc{arc_i}"),
                                    Arc(
                                        source=getattr(b_unit, "waste"),
                                        destination=getattr(
                                            getattr(m.fs, waste_mixer), f"inlet{i}"
                                        ),
                                    ),
                                )
                                arc_i += 1

            # add connection for waste mixer to surface discharge -->
            if "surface_discharge" in unit_list:
                setattr(
                    m.fs,
                    (f"arc{arc_i}"),
                    Arc(
                        source=getattr(m.fs, waste_mixer).outlet,
                        destination=getattr(m.fs, sd_name).inlet,
                    ),
                )
                arc_i += 1
        return m, arc_i, mixer_i

    else:
        return m, arc_i, mixer_i
