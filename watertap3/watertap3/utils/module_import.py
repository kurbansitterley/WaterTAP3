
__all__ = ['get_module']
# from watertap3.wt_units

def get_module(module_name):
    
    # if module_name == 'alum_addition':
    #     from watertap3.unit_models.alum_addition as up

    # if module_name == 'agglom_stacking':
    #     from watertap3.unit_models.agglom_stacking as up

    # if module_name == 'ammonia_addition':
    #     from watertap3.unit_models.ammonia_addition as up

    # if module_name == 'anion_exchange':
    #     from watertap3.unit_models.anion_exchange as up

    # if module_name == 'anion_exchange_epa':
    #     from watertap3.unit_models.anion_exchange_epa as up

    # if module_name == 'anti_scalant_addition':
    #     from watertap3.unit_models.anti_scalant_addition as up

    # if module_name == 'backwash_solids_handling':
    #     from watertap3.unit_models.backwash_solids_handling as up

    # if module_name == 'band_screen':
    #     from watertap3.unit_models.band_screen as up

    if module_name == 'basic_unit':
        from watertap3.unit_models.basic_unit import BasicUnit as up

    # if module_name == 'brine_concentrator':
    #     from watertap3.unit_models.brine_concentrator as up

    if module_name == 'cartridge_filtration':
        from watertap3.unit_models.cartridge_filtration import CartridgeFiltration as up

    # if module_name == 'cation_exchange':
    #     from watertap3.unit_models.cation_exchange as up

    # if module_name == 'caustic_soda_addition':
    #     from watertap3.unit_models.caustic_soda_addition as up

    # if module_name == 'chemical_addition':
    #     from watertap3.unit_models.chemical_addition as up

    if module_name == 'chlorination':
        from watertap3.unit_models.chlorination import Chlorination as up

    # if module_name == 'clarifier':
    #     from watertap3.unit_models.clarifier as up

    # if module_name == 'clearwell':
    #     from watertap3.unit_models.clearwell as up

    # if module_name == 'co2_addition':
    #     from watertap3.unit_models.co2_addition as up

    # if module_name == 'coag_and_floc':
    #     from watertap3.unit_models.coag_and_floc as up

    # if module_name == 'coagulant_addition':
    #     from watertap3.unit_models.coagulant_addition as up

    # if module_name == 'cooling_tower':
    #     from watertap3.unit_models.cooling_tower as up

    # if module_name == 'crystallizer':
    #     from watertap3.unit_models.crystallizer as up

    # if module_name == 'deep_well_injection':
    #     from watertap3.unit_models.deep_well_injection as up

    # if module_name == 'drum_screen':
    #     from watertap3.unit_models.drum_screen as up

    # if module_name == 'electrodialysis_reversal':
    #     from watertap3.unit_models.electrodialysis_reversal as up

    # if module_name == 'evaporation_pond':
    #     from watertap3.unit_models.evaporation_pond as up

    # if module_name == 'ferric_chloride_addition':
    #     from watertap3.unit_models.ferric_chloride_addition as up

    # if module_name == 'filter_press':
    #     from watertap3.unit_models.filter_press as up
    
    # if module_name == 'fixed_bed_gravity_basin':
    #     from watertap3.unit_models.fixed_bed_gravity_basin as up

    # if module_name == 'fixed_bed_pressure_vessel':
    #     from watertap3.unit_models.fixed_bed_pressure_vessel as up

    # if module_name == 'flocculator':
    #     from watertap3.unit_models.flocculator as up

    # if module_name == 'fluidized_bed':
    #     from watertap3.unit_models.fluidized_bed as up
    
    # if module_name == 'gac_gravity':
    #     from watertap3.unit_models.gac_gravity as up

    # if module_name == 'gac_pressure_vessel':
    #     from watertap3.unit_models.gac_pressure_vessel as up
    
    # if module_name == 'gravity_thickener':
    #     from watertap3.unit_models.gravity_thickener as up

    # if module_name == 'grit_chamber':
    #     from watertap3.unit_models.grit_chamber as up

    # if module_name == 'heap_leaching':
    #     from watertap3.unit_models.heap_leaching as up

    # if module_name == 'hydrochloric_acid_addition':
    #     from watertap3.unit_models.hydrochloric_acid_addition as up

    # if module_name == 'ion_exchange':
    #     from watertap3.unit_models.ion_exchange as up

    if module_name == 'iron_and_manganese_removal':
        from watertap3.unit_models.iron_and_manganese_removal import IronAndManganeseRemoval as up

    # if module_name == 'kmno4_addition':
    #     from watertap3.unit_models.kmno4_addition as up

    # if module_name == 'landfill_zld':
    #     from watertap3.unit_models.landfill_zld as up

    # if module_name == 'landfill':
    #     from watertap3.unit_models.landfill as up

    # if module_name == 'lime_addition':
    #     from watertap3.unit_models.lime_addition as up

    # if module_name == 'lime_softening':
    #     from watertap3.unit_models.lime_softening as up

    # if module_name == 'media_filtration':
    #     from watertap3.unit_models.media_filtration as up

    # if module_name == 'micro_screen':
    #     from watertap3.unit_models.micro_screen as up

    # if module_name == 'microfiltration':
    #     from watertap3.unit_models.microfiltration as up

    # if module_name == 'multi_stage_bubble_aeration':
    #     from watertap3.unit_models.multi_stage_bubble_aeration as up

    if module_name == 'treated_water_pumping_station':
        from watertap3.unit_models.treated_water_pumping_station import TreatedWaterPumpingStation as up

    # if module_name == 'nanofiltration':
    #     from watertap3.unit_models.nanofiltration as up

    # if module_name == 'ozone_aop':
    #     from watertap3.unit_models.ozone_aop as up

    # if module_name == 'pac_addition':
    #     from watertap3.unit_models.pac_addition as up

    # if module_name == 'packed_tower_aeration':
    #     from watertap3.unit_models.packed_tower_aeration as up

    if module_name == 'passthrough':
        from watertap3.unit_models.passthrough import Passthrough as up

    # if module_name == 'polymer_addition':
    #     from watertap3.unit_models.polymer_addition as up

    # if module_name == 'rapid_mix':
    #     from watertap3.unit_models.rapid_mix as up

    if module_name == 'reverse_osmosis':
        from watertap3.unit_models.reverse_osmosis import ReverseOsmosis as up

    # if module_name == 'sedimentation':
    #     from watertap3.unit_models.sedimentation as up

    # if module_name == 'sodium_bisulfite_addition':
    #     from watertap3.unit_models.sodium_bisulfite_addition as up

    # if module_name == 'solution_distribution_and_recovery_plant':
    #     from watertap3.unit_models.solution_distribution_and_recovery_plant as up

    # if module_name == 'static_mixer':
    #     from watertap3.unit_models.static_mixer as up

    if module_name == 'storage_tank':
        from watertap3.unit_models.storage_tank import StorageTank as up

    # if module_name == 'sulfuric_acid_addition':
    #     from watertap3.unit_models.sulfuric_acid_addition as up

    # if module_name == 'surface_discharge':
    #     from watertap3.unit_models.surface_discharge as up

    # if module_name == 'sw_onshore_intake':
    #     from watertap3.unit_models.sw_onshore_intake as up

    # if module_name == 'tri_media_filtration':
    #     from watertap3.unit_models.tri_media_filtration as up

    # if module_name == 'ultrafiltration':
    #     from watertap3.unit_models.ultrafiltration as up

    # if module_name == 'uv_aop':
    #     from watertap3.unit_models.uv_aop as up

    # if module_name == 'water_pumping_station':
    #     from watertap3.unit_models.water_pumping_station as up
        
    if module_name == 'well_field':
        from watertap3.unit_models.well_field import WellField as up
    
    # if module_name == 'wire_screen':
    #     from watertap3.unit_models.wire_screen as up

    return up

