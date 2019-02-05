from todloop import TODLoop
from todloop.tod import TODLoader

from cuts import CutSources, CutPlanets, CutPartial, FindJumps, RemoveSyncPickup
from tod import TransformTOD, FouriorTransform, GetDetectors, CalibrateTOD
from analysis import AnalyzeScan, AnalyzeDarkLF

# initialize the pipeline
loop = TODLoop()

# specify the list of tods to go through
loop.add_tod_list('/home/lmaurin/TODLists/2016_ar3_season_nohwp.txt')

# add routines to the pipeline

# add a routine to load tod
loader_params = {
    'output_key': 'tod',
    'repair_pointing': False
}
loop.add_routine(TODLoader(**loader_params))

# add a routine to cut the sources
source_params = {
    'inputs': {
        'tod': 'tod'
    },
    'outputs': {
        'tod': 'tod'
    },
    'tag_source': 'pa3_f90_s16_c10_v1_source',
    'no_noise': True,
    'depot': '/data/actpol/depot'
}
loop.add_routine(CutSources(**source_params))

# add a routine to cut the planets
planets_params = {
    'inputs': {
        'tod': 'tod'
    },
    'outputs': {
        'tod': 'tod',        
    },
    'tag_planet': 'pa3_f90_s16_c10_v1_planet',
    'depot': '/data/actpol/depot',
}
loop.add_routine(CutPlanets(**planets_params))

# add a routine to remove the sync pick up
sync_params = {
    'inputs': {
        'tod': 'tod'
    },
    'outputs': {
        'tod': 'tod',        
    },
    'tag_sync': 'pa3_f90_s16_c10_v1',
    'remove_sync': False,
    'force_sync': False,
    'depot': '/data/actpol/depot',
}
loop.add_routine(RemoveSyncPickup(**sync_params))

# add a routine to cut the glitches
partial_params = {
    'inputs': {
        'tod': 'tod'
    },
    'outputs': {
        'tod': 'tod',        
    },
    'tag_partial': 'pa3_f90_s16_c10_v1_partial',
    'include_mce': True,
    'force_partial': False,
    'glitchp': { 'nSig': 10., 'tGlitch' : 0.007, 'minSeparation': 30, \
                 'maxGlitch': 50000, 'highPassFc': 6.0, 'buffer': 200 },
    'depot': '/data/actpol/depot',
}

loop.add_routine(CutPartial(**partial_params))

# add a routine to transform the TODs
transform_params = {
    'inputs': {
        'tod': 'tod'
    },
    'outputs': {
        'tod': 'tod',        
    },
    'remove_mean': False,
    'remove_median': True,
    'detrend': False,
    'remove_filter_gain': False,
    'n_downsample': 4,  # reduction with 2^4 factor
}
loop.add_routine(TransformTOD(**transform_params))

# add a routine to analyze the scan
scan_params = {
    'inputs': {
        'tod': 'tod'
    },
    'outputs': {
        'scan': 'scan_params',        
    }
}
loop.add_routine(AnalyzeScan(**scan_params))

# add a routine to get the relevant detectors to look at
BASE_DIR = '/data/actpol/actpol_data_shared/ArrayData/2016/ar3/' 
gd_params = {
    'inputs': {
        'tod': 'tod'
    },
    'outputs': {
        'dets': 'dets',        
    },
    'source': 'individual',
    'live': BASE_DIR + 'live_pa3_f90_s16_c10_v4.dict',
    'dark': BASE_DIR + 'dark.dict',
    'exclude': BASE_DIR + 'exclude_pa3_f90_s16_c10_v4.dict'
}
loop.add_routine(GetDetectors(**gd_params))

# add a routine to calibrate DAQ units to pW using flatfield and
# responsivity
cal_params = {
    'inputs': {
        'tod': 'tod',
        'dets': 'dets',
    },
    'outputs': {
        'tod': 'tod',
        'cal': 'calData' 
    },
    'flatfield': "/data/actpol/actpol_data_shared/FlatFields/2015/" + \
                 "ff_actpol3_2015_c9_w1_v2b_mix90-150_it9_actpol3_2015_c9_w2_photon_mix90-150_it7.dict",
    'config': [{
        "type": "depot_cal",
        "depot": "/data/actpol/depot",
        "tag": "pa3_s16_BS",
        "name": "biasstep"
    }, {
        "type": "constant",
        "value": 0.821018,
        "name": "DC bias factor"
    }],
    'forceNoResp': True,
    'calibrateTOD': True,
}
loop.add_routine(CalibrateTOD(**cal_params))

# add a routine to find jumps in TOD
jump_params = {
    'inputs': {
        'tod': 'tod'
    },
    'outputs':{
        'jumps': 'crit_jumps'
    },
    'dsStep': 4,
    'window': 1,
}
loop.add_routine(FindJumps(**jump_params))

# add a routine to perform the fourior transform
fft_params = {
    'inputs': {
        'tod': 'tod'
    },
    'outputs': {
        'tod': 'tod',
        'fft': 'fft_data'
    },            
}
loop.add_routine(FouriorTransform(**fft_params))

# study the dark detectors using LF data
lf_dark_params = {
    'inputs': {
        'tod': 'tod',
        'fft': 'fft_data',
        'dets': 'dets',
        'scan': 'scan_params',
    },
    'outputs': {
        'lf_dark': 'lf_dark',
    },
    'cancelSync': False,
    'doubleMode': False,
    'freqRange': {
        'fmin': 0.017,          
        'fshift': 0.009,
        'band': 0.071,
        'Nwin': 1,
    },
}
loop.add_routine(AnalyzeDarkLF(**lf_dark_params))

# # study the live detectors using LF data
# lf_live_params = {
#     'inputs': {
#         'tod': 'tod',
#         'fft': 'fft_data',
#         'dets': 'dets',
#         'scan': 'scan_params',
#     },
#     'outputs': {
#         'lf_live': 'lf_live',
#     }, 
#     'cancelSync': True,
#     'doubleMode': False,
#     'removeDark': True,
#     'freqRange': {
#         'fmin': 0.017,
#         'fshift': 0.009,
#         'band': 0.071,
#         'Nwin': 10,
#     },
#     'separateFreqs': False,
# }
# loop.add_routine(AnalyzeLiveLF(**lf_live_params))

# run pipeline
loop.run(100,101)
