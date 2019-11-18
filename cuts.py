"""This script aims to implement the cuts in the moby2 framework
without removing the pre-selection. """

import moby2
from todloop import TODLoop

from routines.tod import LoadTOD, CheckTODLength, TransformTOD
from routines.cuts import CutMCE, CutSources, CutPlanets, \
    RemoveSyncPickup, CutPartial, FindPathologies
from routines.misc import FindRebiasTime
from routines.report import PathologyReport

# cut parameter files
cutparam_file = "cutparams_v0.par"
cutParam_file = "cutParams_v0.par"

# load parameters from cutparam file
cutparam = moby2.util.MobyDict.from_file(cutparam_file)
cutParam = moby2.util.MobyDict.from_file(cutParam_file)

# get useful parameters
pathop = cutParam['pathologyParams']
no_noise = not(cutParam.get("fillWithNoise", True))
depot = cutparam['depot']
mask_params = cutParam.get_deep(('source_cuts','mask_params'),{})
shift_params = cutParam.get_deep(('source_cuts', 'mask_shift_generator'))

# initialize loop
loop = TODLoop()

# add list of tod
loop.add_tod_list(cutparam.get('source_scans'))

# find rebias time
config = {
    'config_file': cutparam.get('manifest_conf'),
    'offset': cutparam.get('offset', 0),
    'rebias_wait': cutparam.get('rebias_wait'),
    'IV_wait': cutparam.get('IV_wait')
}
loop.add_routine(FindRebiasTime(**config))

# load tod
loop.add_routine(LoadTOD())

# check whether the tod satisfy the length requirement
config = {
    "fmin": pathop["findPathoParams"]["liveCorrPar"]['freqRange']['fmin'],
    "min_periods": pathop["findPathoParams"].get("minPeriods", 1.)
}
loop.add_routine(CheckTODLength(**config))

# cut mce sources
loop.add_routine(CutMCE(no_noise=no_noise))

# cut sources
config = {
    'inputs': {
        'tod': 'tod'
    },
    'outputs': {
        'tod': 'tod'
    },
    'depot': depot,
    'tag_source': cutparam.get('tag_source'),
    'source_list': cutParam['source_cuts'].get('source_list', None),
    'hdf_source_cuts': cutparam.get('hdf_source_cuts', None),
    'no_noise': no_noise,
    'pointing_par': cutParam['pointing'],
    'mask_params': mask_params,
    'mask_shift_generator': shift_params,
    'write_depot': False,
}
loop.add_routine(CutSources(**config))

# cut planets
config = {
    'inputs': {
        'tod': 'tod'
    },
    'outputs': {
        'tod': 'tod'
    },
    'depot': depot,
    'tag_planet': cutparam.get('tag_planet'),
    'no_noise': no_noise,
    'pointing_par': cutParam['pointing'],
    'mask_params': mask_params,
    'mask_shift_generator': shift_params,
    'write_depot': False,
}
loop.add_routine(CutPlanets(**config))

# remove sync pickup
config = {
    'inputs': {
        'tod': 'tod'
    },
    'outputs': {
        'tod': 'tod'
    },
    'remove_sync': cutparam.get('removeSync', False),
    'force_sync': cutparam.get('forceSync', False),
    'tag_sync': cutparam.get('tag_sync'),
    'depot': depot,
    'write_depot': False
}
loop.add_routine(RemoveSyncPickup(**config))

# partial cuts
config = {
    'inputs': {
        'tod': 'tod'
    },
    'outputs': {
        'tod': 'tod'
    },
    'tag_partial': cutparam.get('tag_partial'),
    'force_partial': cutparam.get('forcePartial'),
    'glitchp': cutParam.get('glitchParams'),
    'include_mce': True,
    'depot': depot,
    'no_noies': no_noise,
    'write': False
}
loop.add_routine(CutPartial(**config))

# transform TOD such as detrend and remove mean
config = {
    'inputs': {
        'tod': 'tod'
    },
    'outputs': {
        'tod': 'tod'
    },
    'remove_mean': False,
    'remove_median': True,
    'detrend': cutparam.get('detrend', False),
    'remove_filter_gain': cutparam.get('remove_filter_gain', True),
    'n_downsample': cutparam.get('n_downsample', 1)
}
loop.add_routine(TransformTOD(**config))

# find pathologies
config = {
    'depot': depot,
    'tag_patho': cutparam.get('tag_patho'),
    'skip_partial': False,
    'force_patho': cutparam.get('forcePatho', False),
    'pathop': pathop
}
loop.add_routine(FindPathologies(**config))

# save pathology report
# FIXME: only support the cutparam file in the same folder
config = {
    'cutparam': cutparam_file,
}
loop.add_routine(PathologyReport(**config))

# run loop
loop.run_parallel(0,4,n_workers=4)
