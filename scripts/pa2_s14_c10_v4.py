from todloop import TODLoop
from todloop.tod import TODLoader

from routines.cuts import CutSources, CutPlanets, CutPartial, FindJumps, RemoveSyncPickup
from routines.tod import TransformTOD, FouriorTransform, GetDetectors, CalibrateTOD
from routines.analysis import AnalyzeScan, AnalyzeDarkLF, AnalyzeLiveLF, GetDriftErrors, \
                              AnalyzeLiveMF, AnalyzeHF
from routines.features import JesseFeatures
from routines.report import Summarize, PrepareDataLabelNew


##############
# parameters #
##############

DEPOT = "/mnt/act3/users/yilun/depot"
actpol_shared = "/mnt/act3/users/yilun/work/actpol_data_shared"
tag = "pa2_s14_c10_v4"

pickle_file = "/mnt/act3/users/yilun/share/pa2/%s_results.pickle" % tag
output_file = "outputs/%s.h5" % tag

n_train = 80
n_validate = 20

#############
# pipeline  #
#############

# initialize the pipelines
train_loop = TODLoop()
validate_loop = TODLoop()
# test_loop = TODLoop()

# specify the list of tods to go through
train_loop.add_tod_list("inputs/%s_train.txt" % tag)
validate_loop.add_tod_list("inputs/%s_validate.txt" % tag)
# test_loop.add_tod_list("inputs/%s_test.txt" % tag)

################################
# add routines to the pipeline #
################################

def add_cut_routines(loop):
    """This function registers a series of common routines for cut
    analysis. This is so that we don't have to keep repeating
    ourselves to register these routines for each data set (train,
    validate, test).
    """

    # add a routine to load tod
    loader_params = {
        'output_key': 'tod',
        'load_opts': {
            'fix_sign': True,
            'repair_pointing': True
        }
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
        'tag_source': '%s_source' % tag,
        'no_noise': True,
        'depot': DEPOT,
        'write_depot': True,
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
        'tag_planet': '%s_planet' % tag,
        'depot': DEPOT,
        'pointing_par': {'source': 'fp_file', \
                         'filename': actpol_shared + "/RelativeOffsets/template_ar2_150201us.txt"
        },
        'mask_params': {
            'radius': (8./60)  #degrees
        },
        'mask_shift_generator': {
            'source':'file',\
            'filename':actpol_shared + '/TODOffsets/tod_offsets_2014_141104_v3.txt',
            'columns': [0,3,4],
            'rescale_degrees': 1./60
        },
        'write_depot': True,
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
        'tag_sync': '%s' % tag,
        'remove_sync': False,
        'force_sync': False,
        'depot': DEPOT,
        'write_depot': True,
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
        'tag_partial': '%s_partial' % tag,
        'include_mce': True,
        'force_partial': False,
        'glitchp': { 'nSig': 10., 'tGlitch' : 0.007, 'minSeparation': 30, \
                     'maxGlitch': 50000, 'highPassFc': 6.0, 'buffer': 200 },
        'depot': DEPOT,
        'write_depot': True,
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
        'n_downsample': 1,  # reduction with 2^n factor
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
    BASE_DIR = '/mnt/act3/users/yilun/work/actpol_data_shared/ArrayData/2015/ar2/'
    gd_params = {
        'inputs': {
            'tod': 'tod'
        },
        'outputs': {
            'dets': 'dets',
        },
        'source': 'individual',
        'live': BASE_DIR + 'live_ext2.dict',
        'dark': BASE_DIR + 'dark_ext.dict',
        'exclude': BASE_DIR + 'exclude.dict'
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
        'flatfield': "/mnt/act3/users/mhasse/shared/actpol_shared_depot/FlatFields/2015/" + \
                     "ff_mixed_ar2_v0_actpol2_2015_c9_v0_photon_it0.dict",
        'config': [{
            "type": "depot_cal",
            "depot": DEPOT,
            "tag": "actpol2_2014_biasstep",
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
            'jumps': 'jumps'
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

    # study the live detectors using LF data
    lf_live_params = {
        'inputs': {
            'tod': 'tod',
            'fft': 'fft_data',
            'dets': 'dets',
            'scan': 'scan_params',
            'dark': 'lf_dark',
            'cal': 'calData'
        },
        'outputs': {
            'lf_live': 'lf_live',
        },
        'cancelSync': True,
        'doubleMode': False,
        'removeDark': True,
        'freqRange': {
            'fmin': 0.017,
            'fshift': 0.009,
            'band': 0.071,
            'Nwin': 10,
        },
        'separateFreqs': False,
        'darkModesParams' : {
            'useDarks': True,
            'useSVD': True,
            'Nmodes': 1,
            'useTherm': False
        },
    }
    loop.add_routine(AnalyzeLiveLF(**lf_live_params))

    # get the drift errors
    de_params = {
        'inputs': {
            'tod': 'tod',
            'fft': 'fft_data',
            'dets': 'dets',
            'scan': 'scan_params',
        },
        'outputs': {
            'drift': 'drift',
        },
        'driftFilter': 0.036,
        'nmodes': 3,
    }
    loop.add_routine(GetDriftErrors(**de_params))

    # study the live detectors in mid-freq
    mf_params = {
        'inputs': {
            'tod': 'tod',
            'fft': 'fft_data',
            'dets': 'dets',
            'scan': 'scan_params',
        },
        'outputs': {
            'mf_live': 'mf_live'
        },
        'midFreqFilter': [0.3, 1.0],
        'nmodes': 8,
    }
    loop.add_routine(AnalyzeLiveMF(**mf_params))

    # study the live and dark detectors in HF
    hf_params = {
        'inputs': {
            'tod': 'tod',
            'fft': 'fft_data',
            'dets': 'dets',
            'scan': 'scan_params',
        },
        'outputs': {
            'hf': 'hf'
        },
        'getPartial': False,
        'highFreqFilter': [9.0, 19.0],
        'nLiveModes': 10,
        'nDarkModes': 3,
        'highOrder': True,
    }
    loop.add_routine(AnalyzeHF(**hf_params))

    # add the routine to compute jesse's features
    params = {
        'inputs': {
            'tod': 'tod',
        },
        'outputs': {
            'results': 'jesse_features',
        }
    }
    loop.add_routine(JesseFeatures(**params))

    # summarize the pickle parameters
    summary_params = {
        'inputs': {
            # calculated features to include in the report
            'features': ['lf_live', 'drift', 'mf_live', 'hf', 'jumps',
                         'jesse_features'],
        },
        'outputs': {
            'report': 'report',
        }
    }
    loop.add_routine(Summarize(**summary_params))

    return loop

#########
# train #
#########

# work on training data
train_loop = add_cut_routines(train_loop)

# save report and TOD data into an h5 file for
# future machine learning pipeline
prepare_params = {
    'inputs': {
        'tod': 'tod',
        'report': 'report',
        'dets': 'dets',
        'fft': 'fft_data',
    },
    'pickle_file': pickle_file,
    'output_file': output_file,
    'group': 'train',
}
train_loop.add_routine(PrepareDataLabelNew(**prepare_params))

# run pipeline for training data
train_loop.run(0, n_train)

############
# validate #
############

# work on validation data
validate_loop = add_cut_routines(validate_loop)

# save report and TOD data into an h5 file for
# future machine learning pipeline
prepare_params.update({
    'group': 'validate'
})
validate_loop.add_routine(PrepareDataLabelNew(**prepare_params))

# run pipeline for validation data
validate_loop.run(0, n_validate)

########
# test #
########

# # work on test data
# test_loop = add_cut_routines(test_loop)

# prepare_params.update({
#     'group': 'test'
# })
# test_loop.add_routine(PrepareDataLabelNew(**prepare_params))

# # run the pipeline for testdata
# test_loop.run(0, n_test)

########
# done #
########
