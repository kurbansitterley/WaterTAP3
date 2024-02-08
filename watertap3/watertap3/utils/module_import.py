
__all__ = ['get_module']
# from watertap3 import wt_units

def get_module(module_name):
    
    if module_name == 'alum_addition':
        from watertap3.unit_models import alum_addition as up

    if module_name == 'agglom_stacking':
        from watertap3.unit_models import agglom_stacking as up

    if module_name == 'ammonia_addition':
        from watertap3.unit_models import ammonia_addition as up

    if module_name == 'anion_exchange':
        from watertap3.unit_models import anion_exchange as up

    if module_name == 'anion_exchange_epa':
        from watertap3.unit_models import anion_exchange_epa as up

    if module_name == 'anti_scalant_addition':
        from watertap3.unit_models import anti_scalant_addition as up

    if module_name == 'backwash_solids_handling':
        from watertap3.unit_models import backwash_solids_handling as up

    if module_name == 'band_screen':
        from watertap3.unit_models import band_screen as up

    if module_name == 'basic_unit':
        from watertap3.unit_models import basic_unit as up

    if module_name == 'brine_concentrator':
        from watertap3.unit_models import brine_concentrator as up

    if module_name == 'cartridge_filtration':
        from watertap3.unit_models import cartridge_filtration as up

    if module_name == 'cation_exchange':
        from watertap3.unit_models import cation_exchange as up

    if module_name == 'caustic_soda_addition':
        from watertap3.unit_models import caustic_soda_addition as up

    if module_name == 'chemical_addition':
        from watertap3.unit_models import chemical_addition as up

    if module_name == 'chlorination':
        from watertap3.unit_models import chlorination as up

    if module_name == 'clarifier':
        from watertap3.unit_models import clarifier as up

    if module_name == 'clearwell':
        from watertap3.unit_models import clearwell as up

    if module_name == 'co2_addition':
        from watertap3.unit_models import co2_addition as up

    if module_name == 'coag_and_floc':
        from watertap3.unit_models import coag_and_floc as up

    if module_name == 'coagulant_addition':
        from watertap3.unit_models import coagulant_addition as up

    if module_name == 'cooling_tower':
        from watertap3.unit_models import cooling_tower as up

    if module_name == 'crystallizer':
        from watertap3.unit_models import crystallizer as up

    if module_name == 'deep_well_injection':
        from watertap3.unit_models import deep_well_injection as up

    if module_name == 'drum_screen':
        from watertap3.unit_models import drum_screen as up

    if module_name == 'electrodialysis_reversal':
        from watertap3.unit_models import electrodialysis_reversal as up

    if module_name == 'evaporation_pond':
        from watertap3.unit_models import evaporation_pond as up

    if module_name == 'ferric_chloride_addition':
        from watertap3.unit_models import ferric_chloride_addition as up

    if module_name == 'filter_press':
        from watertap3.unit_models import filter_press as up
    
    if module_name == 'fixed_bed_gravity_basin':
        from watertap3.unit_models import fixed_bed_gravity_basin as up

    if module_name == 'fixed_bed_pressure_vessel':
        from watertap3.unit_models import fixed_bed_pressure_vessel as up

    if module_name == 'flocculator':
        from watertap3.unit_models import flocculator as up

    if module_name == 'fluidized_bed':
        from watertap3.unit_models import fluidized_bed as up
    
    if module_name == 'gac_gravity':
        from watertap3.unit_models import gac_gravity as up

    if module_name == 'gac_pressure_vessel':
        from watertap3.unit_models import gac_pressure_vessel as up
    
    if module_name == 'gravity_thickener':
        from watertap3.unit_models import gravity_thickener as up

    if module_name == 'grit_chamber':
        from watertap3.unit_models import grit_chamber as up

    if module_name == 'heap_leaching':
        from watertap3.unit_models import heap_leaching as up

    if module_name == 'hydrochloric_acid_addition':
        from watertap3.unit_models import hydrochloric_acid_addition as up

    if module_name == 'ion_exchange':
        from watertap3.unit_models import ion_exchange as up

    if module_name == 'iron_and_manganese_removal':
        from watertap3.unit_models import iron_and_manganese_removal as up

    if module_name == 'kmno4_addition':
        from watertap3.unit_models import kmno4_addition as up

    if module_name == 'landfill_zld':
        from watertap3.unit_models import landfill_zld as up

    if module_name == 'landfill':
        from watertap3.unit_models import landfill as up

    if module_name == 'lime_addition':
        from watertap3.unit_models import lime_addition as up

    if module_name == 'lime_softening':
        from watertap3.unit_models import lime_softening as up

    if module_name == 'media_filtration':
        from watertap3.unit_models import media_filtration as up

    if module_name == 'micro_screen':
        from watertap3.unit_models import micro_screen as up

    if module_name == 'microfiltration':
        from watertap3.unit_models import microfiltration as up

    if module_name == 'multi_stage_bubble_aeration':
        from watertap3.unit_models import multi_stage_bubble_aeration as up

    if module_name == 'municipal_drinking':
        from watertap3.unit_models import municipal_drinking as up

    if module_name == 'nanofiltration':
        from watertap3.unit_models import nanofiltration as up

    if module_name == 'ozone_aop':
        from watertap3.unit_models import ozone_aop as up

    if module_name == 'pac_addition':
        from watertap3.unit_models import pac_addition as up

    if module_name == 'packed_tower_aeration':
        from watertap3.unit_models import packed_tower_aeration as up

    if module_name == 'passthrough':
        from watertap3.unit_models import passthrough as up

    if module_name == 'polymer_addition':
        from watertap3.unit_models import polymer_addition as up

    if module_name == 'rapid_mix':
        from watertap3.unit_models import rapid_mix as up

    if module_name == 'reverse_osmosis':
        from watertap3.unit_models import reverse_osmosis as up

    if module_name == 'sedimentation':
        from watertap3.unit_models import sedimentation as up

    if module_name == 'sodium_bisulfite_addition':
        from watertap3.unit_models import sodium_bisulfite_addition as up

    if module_name == 'solution_distribution_and_recovery_plant':
        from watertap3.unit_models import solution_distribution_and_recovery_plant as up

    if module_name == 'static_mixer':
        from watertap3.unit_models import static_mixer as up

    if module_name == 'storage_tank':
        from watertap3.unit_models import storage_tank as up

    if module_name == 'sulfuric_acid_addition':
        from watertap3.unit_models import sulfuric_acid_addition as up

    if module_name == 'surface_discharge':
        from watertap3.unit_models import surface_discharge as up

    if module_name == 'sw_onshore_intake':
        from watertap3.unit_models import sw_onshore_intake as up

    if module_name == 'tri_media_filtration':
        from watertap3.unit_models import tri_media_filtration as up

    if module_name == 'ultrafiltration':
        from watertap3.unit_models import ultrafiltration as up

    if module_name == 'uv_aop':
        from watertap3.unit_models import uv_aop as up

    if module_name == 'water_pumping_station':
        from watertap3.unit_models import water_pumping_station as up
        
    if module_name == 'well_field':
        from watertap3.unit_models import well_field as up
    
    if module_name == 'wire_screen':
        from watertap3.unit_models import wire_screen as up

    return up

