import numpy as np
from numpy import ma
import scipy.stats.mstats as ms
from scipy import stats as stat
from scipy.cluster.vq import kmeans2
import logging

import moby2
from todloop import Routine

from utils import *


class AnalyzeScan(Routine):
    def __init__(self, **params):
        """This routine analyzes the scan pattern. It takes in an
        TOD and produces scan parameters including the scan freq,
        when it pivot, number of scans. 

        Inputs: 
            tod: TOD data
        Outputs:
            scan_params: 
                T: time per chunk
                pivot: index of pivot
                N: number of chunks
                scan_freq: scan frequency
        """
        Routine.__init__(self)
        self.inputs = params.get('inputs', None)
        self.outputs = params.get('outputs', None)
        self._output_key = params.get('output_key', None)
        self._scan_params = params.get('scan_param', {})
        # self.logger.setLevel(logging.DEBUG)

    def execute(self, store):
        # load tod
        tod = store.get(self.inputs.get('tod'))

        sample_time = (tod.ctime[-1] - tod.ctime[0]) / (tod.ctime.shape[0]-1)
        # analyze scan and save result into a dictionary
        self.logger.info("Analyzing scan pattern...")
        scan = self.analyze_scan(
            np.unwrap(tod.az), sample_time,
            **self._scan_params)

        # get scan frequency
        scan_freq = scan["scan_freq"]

        # get downsample level
        ds = tod.info.downsample_level

        # summary of scan parameters
        scan_params = {
            'T': scan["T"] * ds,
            'pivot': scan["pivot"] * ds,
            'N': scan["N"],
            'scan_freq': scan_freq
        }
        
        self.logger.debug(scan_params)
        store.set(self.outputs.get('scan'), scan_params)

    def analyze_scan(self, az, dt=0.002508, N=50, vlim=0.01, qlim=0.01):
        """Find scan parameters and cuts"""

        # Find no motion

        # compute the 1% and 99% quantiles
        lo, hi = ms.mquantiles(az, (qlim, 1 - qlim))

        # compute the scan speed

        # compute the az steps
        daz = np.diff(az)

        # form the scan speed array
        # note that:
        # daz[0] = az_1 - az_0
        # daz[1] = az_2 - az_1
        # 2*daz[0] - daz[1] = 2 az_1 - 2 az_0 - az_2 + az_1
        #                   = 3 az_1 - 2 az_0 - az_2
        # this is just an estimate of the scan speed at t=0
        v_scan = np.r_[2 * daz[0] - daz[1], daz]

        # smooth the scan speed vector with a simple moving average of
        # length N, the last indexing is to ensure that the size of
        # the array is the same as v_scan
        v_smooth = np.convolve(v_scan,
                               np.ones(N) / N)[(N - 1) / 2:-(N - 1) / 2]

        # estimate the speed using the median of the scan speeds
        speed = np.median(abs(v_smooth))

        # identify when the scanwhen the speed is either too fast
        # or two slow. The minimum speed requirement is specified
        # as a fraction of the median speed using ~vlim~
        stop = abs(v_scan) < vlim * speed
        pick = abs(v_smooth) > 2 * speed

        # exit now in case of a stare TOD
        # i doubt that this is ever going to occur
        if all(stop):
            scan = {
                "az_max": hi,
                "az_min": lo,
                "az_speed": speed,
                "scan_freq": 0.0,
                "az_cuts": None,
                "T": len(az),
                "pivot": 0,
                "N": 1
            }
            return scan

        # Find sections with no scanning by identifying
        # when the speed is below threshold and when the
        # scan range is an outlier
        noscan = stop * (az > lo) * (az < hi)

        # Get scan frequency
        # first calculate the fourior transform
        faz = np.fft.rfft(az - az.mean())
        # identify the highest frequency which corresponds to the
        # scan frequency
        fscan = np.where(abs(faz) == abs(faz).max())[0][0] / dt / len(az)

        # Find turnarounds
        az_min = np.median(az[stop * (az < lo)])
        az_max = np.median(az[stop * (az > hi)])
        # i don't understand this part
        td = abs(lo - az_min) + abs(az_max - hi)

        # Find scan period parameters
        T_scan = int(1. / fscan / dt)  # number of samples in scan period

        # this is kind of arbitrary
        T_ex = int(1.2 * T_scan)

        # find non-stopping part of scan
        onescan = az[~noscan][:T_ex]
        az_min0 = np.min(onescan[stop[~noscan][:T_ex] * (onescan < lo)
                                 * (onescan > lo - td)])
        az_max0 = np.max(onescan[stop[~noscan][:T_ex] * (onescan > hi)
                                 * (onescan < hi + td)])
        imin0 = np.where(az == az_min0)[0][0]
        imax0 = np.where(az == az_max0)[0][0]
        pivot = np.min([imin0, imax0])  # Index of first scan minima or maxima
        N_scan = (len(az) - pivot) / T_scan  # Number of complete scan periods

        # Find cuts
        if hi - lo < 1. * np.pi / 180:
            flag = np.ones_like(az, dtype=bool)
        else:
            flag = (pick + stop) * (az > lo) * (az < hi) + (az < lo - td) + (
                az > hi + td)

        c_vect = moby2.tod.cuts.CutsVector.from_mask(flag).get_buffered(100)

        # return scan parameters
        scan = {
            "az_max": az_max,
            "az_min": az_min,
            "az_speed": speed,
            "scan_freq": fscan,
            "az_cuts": c_vect,
            "T": T_scan,
            "pivot": pivot,
            "N": N_scan
        }
        return scan


class AnalyzeTemperature(Routine):
    def __init__(self, **params):
        """This routine will analyze the temperature of the TOD such as
        measure the mean temperature and thermal drift, and suggest a
        thermalCut

        Inputs:
            tod: TOD data
        Outputs:
            thermal: 
                Temp: mean temperature
                dTemp: thermal drift
                temperatureCut: cuts based on the temperature
        """
        Routine.__init__(self)
        self.inputs = params.get('inputs', None)
        self.outputs = params.get('outputs', None)
        self._channel = params.get('channel', None)
        self._T_max = params.get('T_max', False)
        self._dT_max = params.get('dT_max', None)

    def execute(self, store):
        tod = store.get(self.inputs.get('tod'))

        Temp = None
        dTemp = None
        temperatureCut = False
        
        self.logger.info("Analyzing temperature...")
        if self._channel is None or self._T_max is None \
           or self._dT_max is None:
            pass
        else:
            thermometers = []
            for ch in self._channel:
                thermometer = tod.get_hk(ch, fix_gaps=True)
                if len(np.diff(thermometer).nonzero()[0]) > 0:
                    thermometers.append(thermometer)
            if len(thermometers) > 0:
                thermometers = np.array(thermometers)
                
                # Get thermometer statistics
                th_mean = moby2.tod.remove_mean(data=thermometers)
                th_trend = moby2.tod.detrend_tod(data=thermometers)
                Temp = th_mean[0]
                dTemp = th_trend[1][0] - th_trend[0][0]
                if (Temp > self._T_max) or (abs(dTemp) > self._dT_max):
                    temperatureCut = True
                    
        thermal_results = {
            'Temp': Temp,
            'dTemp': dTemp,
            'temperatureCut': temperatureCut
        }
        
        self.logger.debug(thermal_results)
        store.set(self.outputs.get('thermal'), thermal_results)


class AnalyzeDarkLF(Routine):
    def __init__(self, **params):
        """Analyze the low frequency data and look at the 
        dark detectors only. It calculates properties such
        as corrDark, normDark, gainDark

        Inputs:
            tod: TOD data
            fft: fourior transformed TOD data
            dets: detector lists
            scan: scan parameters
        Outputs:
            results: 
                corrDark: corr for dark detectors
                normDark: norm for dark detectors
                gainDark: gain for dark detectors
        """
        Routine.__init__(self)
        self.inputs = params.get('inputs', None)
        self.outputs = params.get('outputs', None)
        self._dets = params.get('dets', None)
        self._fft_data = params.get('fft_data', None)
        self._tod = params.get('tod', None)
        self._scan = params.get('scan', None)
        self._freqRange = params.get('freqRange', None)
        self._double_mode = params.get('doubleMode', False)
        self._params = params

    def execute(self, store):
        # retrieved relevant data from data store
        tod = store.get(self.inputs.get('tod'))

        fft_data = store.get(self.inputs.get('fft'))
        fdata = fft_data['fdata']
        df = fft_data['df']
        
        sel = store.get(self.inputs.get('dets'))['dark_final']
        scan_freq = store.get(self.inputs.get('scan'))['scan_freq']

        # get the frequency band parameters
        frange = self._freqRange
        fmin = frange.get("fmin", 0.017)
        fshift = frange.get("fshift", 0.009)
        band = frange.get("band", 0.070)
        Nwin = frange.get("Nwin", 1)

        corr = []
        gain = []
        norm = []

        minFreqElem = 16
        # loop over freq windows
        for i in xrange(Nwin):
            # find upper / lower bounds' corresponding index in freq
            # lower bound: fmin + [ fshifts ]
            # upper bound: fmin + band  + [ fshifts ]
            n_l = int(round((fmin + i*fshift)/df))
            n_h = int(round((fmin + i*fshift + band)/df))

            # if there are too few elements then add a few more to
            # have exactly the minimum required
            if (n_h - n_l) < minFreqElem:
                n_h = n_l + minFreqElem

            # perform low frequency analysis
            r = self.lowFreqAnal(fdata, sel, [n_l, n_h], df,
                                 tod.nsamps, scan_freq)

            # append the results to the relevant lists
            corr.append(r["corr"])
            gain.append(np.abs(r["gain"]))
            norm.append(r["norm"])

        # normalize gain
        for g in gain:
            g /= np.mean(g[sel])
            
        gain = np.array(gain)
        
        # give a default gain of 0 for invalid data
        gain[np.isnan(gain)] = 0.

        # use mean as representative values for gain
        mgain = ma.MaskedArray(gain, ~np.array(sel))
        mgain_mean = mgain.mean(axis=0)

        # use max as representative values for corr
        mcorr = ma.MaskedArray(corr, ~np.array(sel))
        mcorr_max = mcorr.max(axis=0)

        # use mean as representative values for norm
        mnorm = ma.MaskedArray(norm, ~np.array(sel))
        mnorm_mean = mnorm.mean(axis=0)
        
        # export the values
        results = {}
        
        results["corrDark"] = mcorr_max.data,
        results["gainDark"] = mgain_mean.data
        results["normDark"] = mnorm_mean.data

        # save to the data store
        store.set(self.outputs.get('results'), results)
        
    def lowFreqAnal(self, fdata, sel, frange, df, nsamps, scan_freq):
        """Find correlations and gains to the main common mode over a
        frequency range
        """
        self.logger.info("Analyzing freqs %s" % frange)
        # get relevant low freq data in the detectors selected
        lf_data = fdata[sel, frange[0]:frange[1]]
        ndet = len(sel)
        res = {}

        # Scan frequency rejection
        if self._params.get("cancelSync", False) and (scan_freq/df > 7):
            i_harm = get_iharm(frange, df, scan_freq,
                               wide=self._params.get("wide",True))
            lf_data[:, i_harm] = 0.0

        # Get correlation matrix
        c = np.dot(lf_data, lf_data.T.conjugate())
        a = np.linalg.norm(lf_data, axis=1)
        aa = np.outer(a,a)
        aa[aa==0.] = 1.
        cc = c/aa

        # Get Norm
        norm = np.zeros(ndet,dtype=float)
        fnorm = np.sqrt(np.abs(np.diag(c)))
        norm[sel] = fnorm*np.sqrt(2./nsamps)
        nnorm = norm/np.sqrt(nsamps)

        # Get Correlations
        u, s, v = np.linalg.svd(lf_data, full_matrices=False)
        corr = np.zeros(ndet)
        if self._double_mode:
            corr[sel] = np.sqrt(abs(u[:,0]*s[0])**2+abs(u[:,1]*s[1])**2)/fnorm
        else:
            corr[sel] = np.abs(u[:,0])*s[0]/fnorm

        # Get Gains
        # data = CM * gain
        gain = np.zeros(ndet, dtype=complex)
        gain[sel] = np.abs(u[:, 0])
        res.update({
            "corr": corr,
            "gain": gain,
            "norm": norm,
            "cc": cc
        })
        return res


class AnalyzeLiveLF(Routine):
    def __init__(self, **params):
        """Analyze the low frequency data and look at the 
        live detectors only. It calculates properties such
        as corrLive, normLive, gainLive and darkRatioLive.

        Inputs:
            tod: TOD data
            fft: fourior transformed TOD data
            dets: detector lists
            scan: scan parameters
            cal: calibration data
        Outputs:
            results: 
                corrLive: corr for live detectors
                normLive: norm for live detectors
                gainLive: gain for live detectors
                darkRatioLive: dark ratio of live detectors
        """
        Routine.__init__(self)
        self.inputs = params.get('inputs', None)
        self.outputs = params.get('outputs', None)
        self._freqRange = params.get('freqRange', None)
        self._separateFreqs = params.get('separateFreqs', False)
        self._full = params.get('fullReport', False)
        self._removeDark = params.get('removeDark', False)
        self._darkModesParams = params.get('darkModesParams', {})
        self._forceResp = params.get("forceResp", True)
        self._params = params

    def execute(self, store):
        # similar to the dark analysis, for more comments please refer to
        # the AnalyzeDarkLF
        
        # retrieve relevant data from store
        # tod, number of detectors and number of samples
        tod = store.get(self.inputs.get('tod'))
        ndets = len(tod.info.det_uid)
        nsamps = tod.nsamps

        # retrieve detector lists
        live = store.get(self.inputs.get('dets'))['live_final']
        dark = store.get(self.inputs.get('dets'))['dark_final']

        # retrieve fft related data
        fft_data = store.get(self.inputs.get('fft'))
        fdata = fft_data['fdata']
        df = fft_data['df']
        nf = fft_data['nf']

        # retrieve scan frequency
        scan_freq = store.get(self.inputs.get('scan'))['scan_freq']

        # retrieve calibration data
        calData = store.get(self.inputs.get('cal'))
        respSel = calData['respSel']
        ff = calData['ff']
        flatfield = calData['flatfield_object']

        # empty list to store the detectors for each frequency bands if
        # that's what we want, or otherwise we will store all detectors here
        fbandSel = []
        fbands = []
        # if we want to treat different frequencies separately
        if self._separateFreqs:
            # gather the different frequency bands
            # i.e. 90GHz, 150GHz, etc
            fbs = np.array(list(set(tod.info.array_data["nom_freq"])))
            fbs = fbs[fbs != 0]
            for fb in fbs:
                # store the live detectors of each frequencies into the respective list
                self.fbandSel.append((tod.info.array_data["nom_freq"] == fb)*live)
                self.fbands.append(str(int(fb)))
        else:
            fbandSel.append(live)
            fbands.append("all")

        # initialize the preselection for live detectors
        preLiveSel = np.zeros(ndets, dtype=bool)
        liveSel = np.zeros(ndets, dtype=bool)

        # initialize vectors to store the statistics for live data
        crit = {}
        crit["darkRatioLive"] = np.zeros(ndets, dtype=float)
        crit["corrLive"] = np.zeros(ndets, dtype=float)
        crit["gainLive"] = np.zeros(ndets, dtype=float)
        crit["normLive"] = np.zeros(ndets, dtype=float)

        # if resp will be used
        if not self._forceResp:
            respSel = None

        # get the frequency band parameters
        frange = self._freqRange
        fmin = frange.get("fmin", 0.017)
        fshift = frange.get("fshift", 0.009)
        band = frange.get("band", 0.070)
        Nwin = frange.get("Nwin", 1)

        # loop over frequency band
        for fbSel,fbn in zip(fbandSel, fbands):
            all_data = []

            sel = []            
            corr = []
            gain = []
            norm = []
            darkRatio = []

            fcm = []
            cm = []
            cmdt = []
            
            minFreqElem = 16
            for i in xrange(Nwin):
                n_l = int(round((fmin + i*fshift)/df))
                n_h = int(round((fmin + i*fshift + band)/df))
                                
                if (n_h - n_l) < minFreqElem:
                    n_h = n_l + minFreqElem

                if self._removeDark:
                    if dark is None:
                        print "ERROR: no dark selection supplied"
                        return 0

                    fcmi, cmi, cmdti = self.getDarkModes(fdata, dark, [n_l,n_h],
                                                         df, nf, nsamps)
                    fcm.append(fcmi)
                    cm.append(cmi)
                    cmdt.append(cmdti)

                r = self.lowFreqAnal(fdata, live, [n_l,n_h], df, nsamps, scan_freq,
                                     fcmodes=fcmi, respSel=respSel, flatfield=flatfield)

                sel.append(live)                
                corr.append(r["corr"])
                gain.append(np.abs(r["gain"]))
                norm.append(r["norm"])
                darkRatio.append(r["ratio"])

                if self._full:
                    all_data.append(r)
                    
            for g in gain:
                g /= np.mean(g[live])
                
            gain = np.array(gain)
            gain[np.isnan(gain)] = 0.


            mgain = ma.MaskedArray(gain,~np.array(sel))
            mgain_mean = mgain.mean(axis=0)
            mcorr = ma.MaskedArray(corr,~np.array(sel))

            mcorr_max = mcorr.max(axis=0)
            mnorm = ma.MaskedArray(norm,~np.array(sel))
            mnorm_mean = mnorm.mean(axis=0)

            # summarize the results so far
            results = {
                "corr": mcorr_max.data,
                "gain": mgain_mean.data,
                "norm": mnorm_mean.data,
            }

            if self._removeDark:
                mdarkRatio = ma.MaskedArray(darkRatio,~np.array(sel))
                mdarkRatio_mean = mdarkRatio.mean(axis=0)
                results['darkRatio'] = mdarkRatio_mean.data

            # update the crit dictionary to output
            crit['corrLive'][fbSel] = results["corr"][fbSel]
            crit['gainLive'][fbSel] = results["gain"][fbSel]
            crit['normLive'][fbSel] = results["norm"][fbSel]

            if results.has_key('darkRatio'):
                crit["darkRatioLive"][fbSel] = results["darkRatio"][fbSel]

        # Undo flatfield correction
        crit["gainLive"] /= np.abs(ff)

        store.set(self.outputs.get('results'), crit)


    def lowFreqAnal(self, fdata, sel, frange, df, nsamps, scan_freq,
                    fcmodes=None, respSel=None, flatfield=None):
        """Find correlations and gains to the main common mode over a
        frequency range
        """
        self.logger.info("Analyzing freqs %s" % frange)
        # get relevant low freq data in the detectors selected
        lf_data = fdata[sel, frange[0]:frange[1]]
        ndet = len(sel)
        res = {}

        # Deproject correlated modes
        if fcmodes is not None:
            data_norm = np.linalg.norm(lf_data,axis=1)
            dark_coeff = []

            # actually do the deprojection here
            for m in fcmodes:
                coeff = np.dot(lf_data.conj(),m)
                lf_data -= np.outer(coeff.conj(),m)
                dark_coeff.append(coeff)

            # Reformat dark coefficients
            if len(dark_coeff) > 0:
                dcoeff = np.zeros([len(dark_coeff),ndet],dtype=complex)
                dcoeff[:,sel] = np.array(dark_coeff)

            # Get Ratio
            ratio = np.zeros(ndet,dtype=float)
            data_norm[data_norm==0.] = 1.
            # after deprojection versus before deprojection
            # this should really be called live ratio than dark ratio
            ratio[sel] = np.linalg.norm(lf_data,axis=1)/data_norm

        # Scan frequency rejection
        if self._params.get("cancelSync",False) and (scan_freq/df > 7):
            i_harm = get_iharm(frange, df, scan_freq,
                               wide=self._params.get("wide",True))
            lf_data[:, i_harm] = 0.0

        # Get correlation matrix
        c = np.dot(lf_data, lf_data.T.conjugate())
        a = np.linalg.norm(lf_data, axis=1)
        aa = np.outer(a,a)
        aa[aa==0.] = 1.
        cc = c/aa

        # Get Norm
        norm = np.zeros(ndet,dtype=float)
        fnorm = np.sqrt(np.abs(np.diag(c)))
        norm[sel] = fnorm*np.sqrt(2./nsamps)
        nnorm = norm/np.sqrt(nsamps)

        # Apply gain ratio in case of multichroic
        if (flatfield is not None) and ("scale" in flatfield.fields):
            scl = flatfield.get_property("scale", det_uid=np.where(sel)[0],
                                         default = 1.)
            lf_data *= np.repeat([scl],lf_data.shape[1],axis=0).T

        # Get Correlations
        u, s, v = np.linalg.svd(lf_data, full_matrices=False )

        corr = np.zeros(ndet)
        if self._params.get("doubleMode", False):
            corr[sel] = np.sqrt(abs(u[:,0]*s[0])**2+abs(u[:,1]*s[1])**2)/fnorm
        else:
            corr[sel] = np.abs(u[:,0])*s[0]/fnorm

        # Get Gains
        # data = CM * gain
        gain = np.zeros(ndet, dtype=complex)
        gain[sel] = np.abs(u[:, 0])
        
        res.update({"corr": corr, "gain": gain, "norm": norm, "dcoeff": dcoeff,
                    "ratio": ratio, "cc": cc})
        
        return res

    def getDarkModes(self, fdata, darkSel, frange, df, nf, nsamps):
        """
        @brief Get dark or thermal modes from dark detectors and thermometer
               data.
        @return correlated modes in frequency and time domain, plus thermometer
                info.
        """
        self.logger.info("Finding dark modes in freqs %s" % frange)        
        n_l, n_h=frange
        fc_inputs = []

        # Dark detector drift
        if self._darkModesParams.get("useDarks", False):
            dark_signal = fdata[darkSel,n_l:n_h].copy()
            fc_inputs.extend(list(dark_signal))

        fc_inputs = np.array(fc_inputs)

        # Normalize modes
        fc_inputs /= np.linalg.norm(fc_inputs, axis=1)[:, np.newaxis]

        # Obtain main svd modes to deproject from data
        if self._darkModesParams.get("useSVD", False):
            Nmodes = self._darkModesParams.get("Nmodes", None)
            u, s, v = np.linalg.svd( fc_inputs, full_matrices=False )
            if Nmodes is None:
                # drop the bottom 10%
                fcmodes = v[s > s.max()/10]
            else:
                fcmodes = v[:Nmodes]
        else:
            fcmodes = fc_inputs

        # Get modes in time domain
        cmodes, cmodes_dt = get_time_domain_modes(
            fcmodes, n_l, nsamps, df)
        cmodes /= np.linalg.norm(cmodes, axis=1)[:, np.newaxis]
        return fcmodes, cmodes, cmodes_dt
        

class GetDriftErrors(Routine):
    def __init__(self, **params):
        """This routine obtains the pickle parameter DELive by performing a
        high frequency analysis on the slow modes
        
        Inputs:
            tod: TOD data
            dets: detector lists
            fft: fft data
            scan: scan parameters
        Outputs:
            results:
                DELive: drift error for live detectors
        """
        Routine.__init__(self)
        self.inputs = params.get('inputs', None)
        self.outputs = params.get('outputs', None)
        self._driftFilter = params.get('driftFilter', None)
        self._nmodes = params.get('nmodes', 1)

    def execute(self, store):
        tod = store.get(self.inputs.get('tod'))
        nsamps = tod.nsamps

        live = store.get(self.inputs.get('dets'))['live_final']

        fft_data = store.get(self.inputs.get('fft'))
        fdata = fft_data['fdata']
        df = fft_data['df']

        scan_freq = store.get(self.inputs.get('scan'))['scan_freq']
        nmodes = self._nmodes

        # find the range of frequencies of interests
        n_l = 1
        n_h = nextregular(int(round(self._driftFilter/df))) + 1

        # get drift errors
        ndets = len(live)
        hf_data = fdata[live, n_l:n_h]

        # remove first [nmodes] common modes
        if nmodes > 0:
            # find the correlation between different detectors
            c = np.dot(hf_data, hf_data.T.conjugate())

            # find the first few common modes in the detectors and
            # deproject them
            u, w, v = np.linalg.svd(c, full_matrices = 0)
            kernel = v[:nmodes]/np.repeat([np.sqrt(w[:nmodes])],len(c),axis=0).T
            modes = np.dot(kernel, hf_data)
            coeff = np.dot(modes, hf_data.T.conj())
            hf_data -= np.dot(coeff.T.conj(), modes)

        # compute the rms for the detectors
        rms = np.zeros(ndets)
        rms[live] = np.sqrt(np.sum(abs(hf_data)**2,axis=1)/hf_data.shape[1]/nsamps)

        results = {
            "DELive": rms
        }

        store.set(self.outputs.get('results'), results)


class AnalyzeLiveMF(Routine):
    def __init__(self, **params):
        """This routine looks at the mid-frequency and perform a 
        high-freq like analysis to get the pickle parameter MFE

        Inputs:
            tod: TOD data
            dets: detector lists
            fft: fft data
            scan: scan parameters
        Outputs:
            results:
                MFELive: mid frequency error for live detectors
        """
        Routine.__init__(self)
        self.inputs = params.get('inputs', None)
        self.outputs = params.get('outputs', None)
        self._midFreqFilter = params.get("midFreqFilter", None)
        self._nmodes = params.get("nmodes", 1)

    def execute(self, store):
        tod = store.get(self.inputs.get('tod'))
        nsamps = tod.nsamps

        live = store.get(self.inputs.get('dets'))['live_final']

        fft_data = store.get(self.inputs.get('fft'))
        fdata = fft_data['fdata']
        df = fft_data['df']

        scan_freq = store.get(self.inputs.get('scan'))['scan_freq']
        nmodes = self._nmodes

        # get the frequency range to work on
        n_l = int(self._midFreqFilter[0]/df)
        n_h = int(self._midFreqFilter[1]/df)

        # get drift errors
        ndets = len(live)
        hf_data = fdata[live, n_l:n_h]
        
        # remove first [nmodes] common modes
        if nmodes > 0:
            self.logger.info("Deprojecting %d modes" % nmodes)
            # find the correlation between different detectors
            c = np.dot(hf_data, hf_data.T.conjugate())

            # find the first few common modes in the detectors and
            # deproject them
            u, w, v = np.linalg.svd(c, full_matrices = 0)
            kernel = v[:nmodes]/np.repeat([np.sqrt(w[:nmodes])],len(c),axis=0).T
            modes = np.dot(kernel, hf_data)
            coeff = np.dot(modes, hf_data.T.conj())
            hf_data -= np.dot(coeff.T.conj(), modes)

        # compute the rms for the detectors
        rms = np.zeros(ndets)
        rms[live] = np.sqrt(np.sum(abs(hf_data)**2,axis=1)/hf_data.shape[1]/nsamps)

        results = {
            "MFELive": rms
        }

        store.set(self.outputs.get('results'), results)


class AnalyzeHF(Routine):
    def __init__(self, **params):
        """This routine analyzes both live and dark detectors in
        the high frequency band. It calculates properties such as
        rmsLive, skewLive, kurtLive and their dark equivalent

        Inputs:
            tod: TOD data
            dets: detector lists
            fft: fft data
            scan: scan parameters
        Outputs:
            results:
                rmsLive, skewLive, kurtLive, rmsDark, skewLive, 
                kurtLive
        """
        Routine.__init__(self)
        self.inputs = params.get('inputs', None)
        self.outputs = params.get('outputs', None)
        self._getPartial = params.get('getPartial', False)
        self._highFreqFilter = params.get('highFreqFilter', None)
        self._nmodes_live = params.get('nLiveModes', 1)
        self._nmodes_dark = params.get('nDarkModes', 1)
        self._highOrder = params.get('highOrder', False)
        self._params = params

    def execute(self, store):
        # load relevant data from the data store
        tod = store.get(self.inputs.get('tod'))
        nsamps = tod.nsamps
        ndets = len(tod.info.det_uid)

        live = store.get(self.inputs.get('dets'))['live_final']
        dark = store.get(self.inputs.get('dets'))['dark_final']

        fft_data = store.get(self.inputs.get('fft'))
        fdata = fft_data['fdata']
        df = fft_data['df']

        scan = store.get(self.inputs.get('scan'))
        nmodes_live = self._nmodes_live
        nmodes_dark = self._nmodes_dark

        # get the range of frequencies to work with
        n_l = int(round(self._highFreqFilter[0]/df))
        n_h = int(round(self._highFreqFilter[1]/df))

        # make sure that n_h is a number that's optimized in fft
        n_h = nextregular(n_h-n_l) + n_l

        # empty dictionary to store the results
        results = {}

        # if partial is labeled the analysis will be carried out
        # for each individual scan between the turnning points
        # TODO: There is still some problem with partial analysis
        if not(self._getPartial):
            self.logger.info("Performing non-partial analysis")
            rms, skewt, kurtt = self.highFreqAnal(fdata, live, [n_l,n_h], nsamps,
                                                  highOrder=self._highOrder,
                                                  nmodes=nmodes_live)
        else:
            self.logger.info("Performing partial analysis")
            rms, skewt, kurtt, prms, pskewt, pkurtt = self.highFreqAnal(fdata, live, 
                                                                        [n_l,n_h],
                                                                        nsamps,
                                                                        nmodes=nmodes_live,
                                                                        highOrder=self._highOrder,
                                                                        scanParams=scan)

            # store the statistics for partial
            results["partialRMSLive"] = np.zeros([ndets, scan["N"]])
            results["partialSKEWLive"] = np.zeros([ndets, scan["N"]]) 
            results["partialKURTLive"] = np.zeros([ndets, scan["N"]]) 
            results["partialSKEWPLive"] = np.zeros([ndets, scan["N"]])
            results["partialKURTPLive"] = np.zeros([ndets, scan["N"]])

            results["partialRMSLive"][live] =  prms
            results["partialSKEWLive"][live] = pskewt.T[:, 0]
            results["partialKURTLive"][live] = pkurtt.T[:, 0]
            results["partialSKEWPLive"][live] = pskewt.T[:, 1]
            results["partialKURTPLive"][live] = pkurtt.T[:, 1]
            
        # store the statistics for global
        results["rmsLive"] = rms
        results["skewLive"] = np.zeros(ndets)
        results["kurtLive"] = np.zeros(ndets)
        results["skewpLive"] = np.zeros(ndets)
        results["kurtpLive"] = np.zeros(ndets)
        results["skewLive"][live] = skewt[0]
        results["kurtLive"][live] = kurtt[0]
        results["skewpLive"][live] = skewt[1]
        results["kurtpLive"][live] = kurtt[1]

        # analyze the dark detectors for the same frequency range
        self.logger.info("Analyzing dark detectors")
        rms = self.highFreqAnal(fdata, dark, [n_l,n_h], nsamps, nmodes=nmodes_dark, 
                                highOrder=False)

        results["rmsDark"] = rms

        store.set(self.outputs.get('results'), results)

    def highFreqAnal(self, fdata, sel, frange, nsamps, nmodes=0, highOrder=False,
                     scanParams=None):
        """
        @brief Find noise RMS, skewness and kurtosis over a frequency band
        """
        self.logger.info("Analyzing freqs %s" % frange)
        ndet = len(sel)

        # get the high frequency fourior modes
        hf_data = fdata[sel, frange[0]:frange[1]]

        if nmodes > 0:
            self.logger.info("Deprojecting %d modes" % nmodes)
            # find the correlation between different detectors
            c = np.dot(hf_data,hf_data.T.conjugate())

            # find the first few common modes in the detectors and
            # deproject them
            u, w, v = np.linalg.svd(c, full_matrices = 0)
            kernel = v[:nmodes]/np.repeat([np.sqrt(w[:nmodes])],len(c),axis=0).T
            modes = np.dot(kernel,hf_data)
            coeff = np.dot(modes,hf_data.T.conj())
            hf_data -= np.dot(coeff.T.conj(),modes)

        # compute the rms for the detectors
        rms = np.zeros(ndet)
        rms[sel] = np.sqrt(np.sum(abs(hf_data)**2,axis=1)/hf_data.shape[1]/nsamps)

        # if we are interested in high order effects, skew and kurtosis will
        # be calculated here
        if highOrder:
            hfd, _ = get_time_domain_modes(hf_data, 1, nsamps)
            skewt = stat.skewtest(hfd,axis=1)
            kurtt = stat.kurtosistest(hfd,axis=1)
            if scanParams is not None:
                T = scanParams["T"]
                pivot = scanParams["pivot"]
                N = scanParams["N"]
                f = float(hfd.shape[1])/nsamps
                t = int(T*f); p = int(pivot*f)
                prms = []; pskewt = []; pkurtt = []
                for c in xrange(N):
                    # i see this as calculating the statistics for each
                    # swing (no turning part) so the statistics is not
                    # affected by the scan
                    prms.append(hfd[:,c*t+p:(c+1)*t+p].std(axis=1))
                    pskewt.append(stat.skewtest(hfd[:,c*t+p:(c+1)*t+p],axis=1)) 
                    pkurtt.append(stat.kurtosistest(hfd[:,c*t+p:(c+1)*t+p],axis=1)) 
                prms = np.array(prms).T
                pskewt = np.array(pskewt)
                pkurtt = np.array(pkurtt)
                return (rms, skewt, kurtt, prms, pskewt, pkurtt)
            else:
                return (rms, skewt, kurtt)
        else:
            return rms        
