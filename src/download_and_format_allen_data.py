import os
from pathlib import Path
import loading_classes as cls

directory = os.path.join(os.getcwd(), 'src')

# read in object
session_id = [799864342,798911424,797828357,791319847,773418906,763673393,762602078,762120172,
             761418226,760693773,760345702,759883607,758798717,757970808,757216464,756029989,
             755434585,754829445,754312389,751348571,750749662,750332458,746083955,744228101,
             743475441,742951821,739448407,737581020,732592105,721123822,719161530,715093703]
stimulus = cls.Stimulus()
stimulus.build_region_unit_config_trial_spikebins_multi_parallel(session_id)

# stimulus.save_session()

# # save object
# dirr = os.path.join(directory, 'pickle_files')
# os.makedirs(dirr, exist_ok=True)
# stimulus.write_sessions_to_disk(dirr, 'ALL_SESSIONS')

# # load object
# dirr = os.path.join(directory, 'pickle_files')
# saved_stimulus = cls.load_sessions_from_disk(dirr, 'ALL_SESSIONS')
# stimulus = cls.Stimulus(old_object=saved_stimulus)

# read out trial-level psth data for analysis in R.
curr_sesh = stimulus._current_session
regions = curr_sesh._loaded_regions
orients = curr_sesh._loaded_orientations
freqs = curr_sesh._loaded_frequencies
psth_dir = os.path.join(os.path.expanduser("~"), 'RStudioProjects', 'rNeuroPixel', 'rDataset', 'units_data_'+str(curr_sesh._id))
os.makedirs(psth_dir, exist_ok=True)
configs = [(x, y) for x in orients for y in freqs]
for orient, freq in configs:
    presentation_ids = curr_sesh._presentations[(curr_sesh._presentations.orientation == orient) &
                                               (curr_sesh._presentations.temporal_frequency == freq)].index.values
    for region in regions:
        units = curr_sesh._units[curr_sesh._units.ecephys_structure_acronym ==
                                 curr_sesh._structure_map[region]].index.values
        save_psth_dir = os.path.join(psth_dir, str(freq)+'-'+str(orient), region)
        os.makedirs(save_psth_dir, exist_ok=True)
        for unit in units:
            unit_data = curr_sesh._data.sel(unit_id=unit, stimulus_presentation_id=presentation_ids).\
                to_pandas().transpose()
            unit_data.to_csv(os.path.join(save_psth_dir, str(unit) + '.csv'), index=None)
