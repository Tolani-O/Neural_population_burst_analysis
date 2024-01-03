#Version with no imports drifting gratings contrasts.

import os
import numpy as np
from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
import pickle
from multiprocessing import Pool
import traceback
import platform


class SessionData:
    _structure_map = {
        'V1': 'VISp',
        'LM': 'VISl',
        'AL': 'VISal',
        'AM': 'VISam',
        'PM': 'VISpm',
        'RL': 'VISrl',
        'TH': 'LGd'
    }
    _id = None
    _session = None
    _presentations = None
    _units = None
    _data = None
    _unit_counts = {}
    _loaded_regions = None
    _loaded_conditions = None
    _time_step = None
    _time_bin = None

    def __init__(self, session, id):
        self._id = id
        self._session = session

    def getdata(self, regions, frequencies, orientations, time_step, time_bin):
        print('STARTED: Loading Session:', self._id, 'with worker:', os.getpid(), 'from parent', os.getppid())
        self._time_step = time_step
        self._time_bin = time_bin
        self._presentations = self._session.get_stimulus_table("drifting_gratings")
        self._presentations = self._presentations[
            self._presentations.orientation.isin(orientations) &
            self._presentations.temporal_frequency.isin(frequencies)]
        conditions = list(self._presentations.groupby(['orientation','temporal_frequency']).groups.keys())
        self._loaded_conditions = [tuple(int(x) for x in y) for y in conditions]
        ccf_regions = [self._structure_map[x] for x in regions]
        self._units = self._session.units[self._session.units.ecephys_structure_acronym.isin(ccf_regions)]
        bin_edges = np.arange(time_bin[0], time_bin[1] + time_step, time_step)
        self._data = self._session.presentationwise_spike_counts(
            stimulus_presentation_ids=self._presentations.index.values,
            bin_edges=bin_edges,
            unit_ids=self._units.index.values)
        unit_cts = self._units.ecephys_structure_acronym.value_counts()
        self._unit_counts = {list(self._structure_map.keys())[list(self._structure_map.values()).index(x)]: unit_cts[x] for x in list(unit_cts.keys())}
        self._loaded_regions = list(self._unit_counts)
        print('COMPLETED: Loading Session:', self._id, 'with worker:', os.getpid(), 'from parent', os.getppid())


class Stimulus:
    _orientations = [0, 45, 90, 135, 180, 225, 270, 315]
    _frequencies = [1, 2, 4, 8, 15]
    _cache = None
    _sessions = dict()
    _current_session = None
    _manifest_path = ''
    _output_folder = ''
    _log_file = ''

    def __init__(self, old_object=None):
        environment = platform.system()
        if environment=='Darwin':
            self._manifest_path = os.path.join(os.path.expanduser('~'), 'PycharmProjects', 'NeuroPixel', 'Dataset', 'ecephys_cache_dir')
            self._output_folder = os.path.join(os.path.expanduser('~'), 'RStudioProjects', 'rNeuroPixel', 'rDataset')
        elif environment=='Linux':
            self._manifest_path = os.path.join(os.path.expanduser('~'), 'data', 'python', 'ecephys_cache_dir')
            self._output_folder = os.path.join(os.path.expanduser('~'), 'data', 'r', 'dataset')
        os.makedirs(self._manifest_path, exist_ok=True)
        os.makedirs(self._output_folder, exist_ok=True)
        self._manifest_path = os.path.join(self._manifest_path, 'manifest.json')
        self._log_file = os.path.join(self._output_folder, 'log.txt')
        if old_object is None:
            self._cache = EcephysProjectCache.from_warehouse(manifest=self._manifest_path)
        else:
            self.hydrate_object(old_object)

    def hydrate_object(self, old_object):
        self._cache = old_object._cache
        self._sessions = old_object._sessions
        self._current_session = old_object._current_session

    def get_session_single(self, session_id):
        print('STARTED: Getting Session:', session_id, 'with worker:', os.getpid(), 'from parent', os.getppid())
        if session_id not in self._sessions:
            session = self._cache.get_session_data(session_id)
            self._sessions[session_id] = SessionData(session, session_id)
        self._current_session = self._sessions[session_id]
        print('COMPLETED: Getting Session:', session_id, 'with worker:', os.getpid(), 'from parent', os.getppid())

    def build_region_unit_config_trial_spikebins_multi(self, session_id_multi):
        for session_id in session_id_multi:
            self.build_region_unit_config_trial_spikebins_single(session_id, disk=False, cache=True, log=False)

    def build_region_unit_config_trial_spikebins_multi_parallel(self, session_id_multi):
        with open(self._log_file, 'w') as f:
            f.write('LOGS\n\n')
        with Pool() as pool:
            pool.map(self.build_region_unit_config_trial_spikebins_single, session_id_multi)

    def build_region_unit_config_trial_spikebins_single(self, session_id, regions=None, frequencies=None,
                                                        orientations=None, time_step=None, time_bin=None,
                                                        disk=False, cache=False, log=False):
        try:
            if session_id is not None:  # None here would need to be explicits passed in
                self.get_session_single(session_id)
            if time_bin is None:
                time_bin = [0, 2]
            if time_step is None:
                time_step = 0.002
            if regions is None:
                regions = list(self._current_session._structure_map.keys())
            if frequencies is None:
                frequencies = self._frequencies
            if orientations is None:
                orientations = self._orientations
            self._current_session.getdata(regions, frequencies, orientations, time_step, time_bin)
            self.save_current_session(disk, cache)
        except:
            msg = 'FAILED: Session: ' + str(self._current_session._id) + ' with worker: ' + str(os.getpid()) + \
                  ' from parent ' + str(os.getppid())
            print(msg)
            traceback.print_exc()
            if log:
                with open(self._log_file, 'a') as f:
                    f.write(msg+'\n')
                    traceback.print_exc(-1, f)
                    f.write('\n\n')

    def save_current_session(self, disk=True, cache=False):
        print('STARTED: Storing Session:', self._current_session._id, 'with worker:', os.getpid(), 'from parent', os.getppid())
        if disk:
            self.write_current_session_to_disk()
        if cache:
            self._sessions[self._current_session._id] = self._current_session
        print('COMPLETED: Storing Session:', self._current_session._id, 'with worker:', os.getpid(), 'from parent', os.getppid())

    def write_current_session_to_disk(self):
        curr_sesh = self._current_session
        regions = curr_sesh._loaded_regions
        conditions = curr_sesh._loaded_conditions

        psth_dir = os.path.join(self._output_folder, 'units_data_'+str(curr_sesh._id))
        os.makedirs(psth_dir, exist_ok=True)

        for orient, freq in conditions:
            presentation_ids = curr_sesh._presentations[(curr_sesh._presentations.orientation == orient) &
                                                        (curr_sesh._presentations.temporal_frequency == freq)].index.values
            for region in regions:
                units = curr_sesh._units[curr_sesh._units.ecephys_structure_acronym ==
                                         curr_sesh._structure_map[region]].index.values
                save_psth_dir = os.path.join(psth_dir, str(freq)+'-'+str(orient), region)
                os.makedirs(save_psth_dir, exist_ok=True)
                for unit in units:
                    unit_data = curr_sesh._data.sel(unit_id=unit, stimulus_presentation_id=presentation_ids). \
                        to_pandas().transpose()
                    unit_data.to_csv(os.path.join(save_psth_dir, str(unit) + '.csv'), index=None)

    def pickle_all_sessions(self, directory, file_name):
        with open(os.path.join(directory, file_name + '.p'), 'wb') as object_file:
            pickle.dump(self, object_file, pickle.HIGHEST_PROTOCOL)


def load_sessions_from_disk(directory, file_name):
    with open(os.path.join(directory, file_name + '.p'), 'rb') as object_file:
        p_object = pickle.load(object_file)
    return p_object
