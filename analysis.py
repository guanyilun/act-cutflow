import numpy as np
from numpy import ma
import scipy.stats.mstats as ms
from scipy.cluster.vq import kmeans2

import moby2
from todloop import Routine

from utils import *


class AnalyzeScan(Routine):
    def __init__(self, **params):
        """This routine analyzes the scan pattern"""
        Routine.__init__(self)
        self._input_key = params.get('input_key', None)
        self._output_key = params.get('output_key', None)
        self._scan_params = params.get('scan_param', {})

    def execute(self, store):
        # load tod
        tod = store.get(self._input_key)

        sample_time = (tod.ctime[-1] - tod.ctime[0]) / (tod.ctime.shape[0]-1)
        # analyze scan and save result into a dictionary
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
        
        self.logger.info(scan_params)
        store.set(self._output_key, scan_params)

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
        """This routine will analyze the temperature of the TOD"""
        Routine.__init__(self)
        self._input_key = params.get('input_key', None)
        self._output_key = params.get('output_key', None)
        self._channel = params.get('channel', None)
        self._T_max = params.get('T_max', False)
        self._dT_max = params.get('dT_max', None)

    def execute(self, store):
        """
        @brief Measure the mean temperature and thermal drift, and
               suggest a thermalCut
        @return mean temperature, thermal drift and thermal cut flag
        """
        tod = store.get(self._input_key)

        Temp = None
        dTemp = None
        temperatureCut = False
        
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
        
        self.logger.info(thermal_results)
        store.set(self._output_key, thermal_results)


class FouriorTransform(Routine):
    def __init__(self, **params):
        Routine.__init__(self)
        self._input_key = params.get('input_key', None)
        self._output_key = params.get('output_key', None)
        self._fft_data = params.get('fft_data', None)

    def execute(self, store):
        tod = store.get(self._input_key)

        # first de-trend tod 
        trend = moby2.tod.detrend_tod(tod)

        # find the next regular, this is to make fft faster
        nf = nextregular(tod.nsamps)
        fdata = np.fft.rfft(tod.data, nf)

        # time and freq units
        dt = (tod.ctime[-1]-tod.ctime[0])/(tod.nsamps-1)
        df = 1./(dt*nf)

        # summarize fft data
        fft_data = {
            'trend': trend,
            'fdata': fdata,
            'dt': dt,
            'df': df,
            'nf': nf
        }

        # store data into data store
        store.set(self._output_key, tod)
        store.set(self._fft_data, fft_data)


class AnalyzeDarkLF(Routine):
    def __init__(self, **params):
        self._dets = params.get('dets', None)
        self._fft_data = params.get('fft_data', None)
        self._tod = params.get('tod', None)
        self._output_key = params.get('output_key', None)
        self._scan = params.get('scan', None)
        self._freqRange = params.get('freqRange', None)
        self._params = params

    def execute(self, store):
        # retrieved relevant data from data store
        tod = store.get(self._tod)
        fft_data = store.get(self._fft_data)
        fdata = fft_data['fdata']
        df = fft_data['df']
        sel = store.get(self._dets)['dark_final']
        scan_freq = store.get(self._scan)['scan_freq']

        # get the frequency band parameters
        frange = self._freqRange
        fmin = frange.get("fmin", 0.017)
        fshift = frange.get("fshift", 0.009)
        band = frange.get("band", 0.070)
        Nwin = frange.get("Nwin", 1)

        psel = []
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
            psel.append(r["preSel"])
            corr.append(r["corr"])
            gain.append(np.abs(r["gain"]))
            norm.append(r["norm"])

        # count the decision for each frequency window, for example
        # if we looked at 10 frequency window, psel may look like
        # [S S U S S U S U S U] with S means selected as good and
        # u means unselected (bad). Here we want to count the number
        # of "votes" saying this detector is good
        psel = np.sum(psel, axis=0)

        # get the highest amount of score that the detectors receive
        # and use this as a threshold to judge the rest of the
        # detectors
        Nmax = psel.max()

        # For any detectors who are "voted" as good more than half
        # of the maximum "votes" will be selected as good
        psel50 = psel >= Nmax/2.

        # Normalize gain by the average gain of a good selection of
        # detectors, here this selection is given by psel50 AND presel
        for g, s in zip(gain, psel):
            g /= np.mean(g[psel50*s])
            
        gain = np.array(gain)
        
        # give a default gain of 0 for invalid data
        gain[np.isnan(gain)] = 0.

        # use mean as representative values for gain
        mgain = ma.MaskedArray(gain, ~np.array(psel))
        mgain_mean = mgain.mean(axis=0)

        # use max as representative values for corr
        mcorr = ma.MaskedArray(corr, ~np.array(psel))
        mcorr_max = mcorr.max(axis=0)

        # use mean as representative values for norm
        mnorm = ma.MaskedArray(norm, ~np.array(psel))
        mnorm_mean = mnorm.mean(axis=0)
        
        # export the values
        results = {}
        
        results["corrDark"] = mcorr_max.data,
        results["gainDark"] = mgain_mean.data
        results["normDark"] = mnorm_mean.data
        results["darkSel"] = psel50.copy()  # not sure why copy is needed
        store.set(self._output_key, results)
        
    def lowFreqAnal(self, fdata, sel, frange, df, nsamps, scan_freq):
        """Find correlations and gains to the main common mode over a
        frequency range
        """
        # get relevant low freq data in the detectors selected
        lf_data = fdata[sel, frange[0]:frange[1]]
        ndet = len(sel)
        res = {}

        # Apply sine^2 taper to data
        if self._params.get("useTaper", False):
            taper = get_sine2_taper(frange, edge_factor = 6)
            lf_data *= np.repeat([taper],len(lf_data),axis=0)

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
        ppar = self._params.get("presel",{})
        norm = np.zeros(ndet,dtype=float)
        fnorm = np.sqrt(np.abs(np.diag(c)))
        norm[sel] = fnorm*np.sqrt(2./nsamps)
        nnorm = norm/np.sqrt(nsamps)

        # get a range of valid norm values 
        nlim = ppar.get("normLimit",[0.,1e15])
        if np.ndim(nlim) == 0:
            nlim = [0, nlim]
        normSel = (nnorm > nlim[0])*(nnorm < nlim[1])
        
        # check which preselection is specified
        presel_method = ppar.get("method", "median")
        if presel_method is "median":
            sl = presel_by_median(cc, sel=normSel[sel], **presel_params)
            res["groups"] = None
            
        elif presel_method is "groups":
            G, ind, ld, smap = group_detectors(cc, sel=normSel[sel], **presel_params)
            sl = np.zeros(cc.shape[1], dtype=bool)
            sl[ld] = True
            res["groups"] = {
                "G": G,
                "ind": ind,
                "ld": ld,
                "smap": smap
            }
        else:
            raise "ERROR: Unknown preselection method"

        # The number of sels are just overwhelmingly confusing
        # To clarify for myself,
        # - normSel: selects the detectors with good norm
        # - sel: the initial selection of detectors specified
        #        for this case it is selection of dark detectors
        # - sl: is the preselected detectors from the median
        #       or group methods
        # Here it's trying to apply the preselection to the
        # dark selection
        preSel = sel.copy()
        preSel[sel] = sl
        
        # Get Correlations
        u, s, v = np.linalg.svd(lf_data[sl], full_matrices=False )
        corr = np.zeros(ndet)
        if par.get("doubleMode", False):
            corr[preSel] = np.sqrt(abs(u[:,0]*s[0])**2+abs(u[:,1]*s[1])**2)/fnorm[sl]
        else:
            corr[preSel] = np.abs(u[:,0])*s[0]/fnorm[sl]

        # Get Gains
        # data = CM * gain
        gain = np.zeros(ndet, dtype=complex)
        gain[preSel] = np.abs(u[:, 0])
        res.update({"preSel": preSel, "corr": corr, "gain": gain, "norm": norm, 
                    "cc": cc, "normSel": normSel})
        
        return res


class AnalyzeLiveLF(Routine):
    def __init__(self, **params):
        Routine.__init__(self)
        self._dets = params.get('dets', None)
        self._fft_data = params.get('fft_data', None)
        self._tod = params.get('tod', None)
        self._output_key = params.get('output_key', None)
        self._scan = params.get('scan', None)
        self._freqRange = params.get('freqRange', None)
        self._separateFreqs = params.get('separateFreqs', False)
        self._dark = params.get('dark_results', None)
        self._full = params.get('fullReport', False)
        self._removeDark = params.get('removeDark', False)
        self._params = params

    def execute(self, store):
        # similar to the dark analysis, for more comments please refer to
        # the AnalyzeDarkLF
        
        # retrieve relevant data from store
        tod = store.get(self._tod)
        live = store.get(self._dets)['live_final']
        ndets = len(tod.info.det_uid)
        fdata = store.get(self._fft_data)['fdata']
        df = store.get(self._fft_data)['df']
        nf = store.get(self._fft_data)['nf']
        scan_freq = store.get(self._scan)['scan_freq']
        darkSel = store.get(self._dark)['darkSel']
        
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
            self.fbandSel.append(live)
            self.fbands.append("all")

        # initialize the preselection for live detectors
        self.preLiveSel = np.zeros(ndets, dtype=bool)
        self.liveSel = np.zeros(ndets, dtype=bool)

        # initialize the vectors to store the statistics or live data
        self.crit["darkRatioLive"] = np.zeros(ndets, dtype=float)
        self.crit["corrLive"] = np.zeros(ndets, dtype=float)
        self.crit["gainLive"] = np.zeros(ndets, dtype=float)
        self.crit["normLive"] = np.zeros(ndets, dtype=float)

        # if resp will be used
        if not(par[parTag].get("forceResp",True)):
            respSel = None
                
        # get the frequency band parameters
        frange = self._freqRange
        fmin = frange.get("fmin", 0.017)
        fshift = frange.get("fshift", 0.009)
        band = frange.get("band", 0.070)
        Nwin = frange.get("Nwin", 1)

        # initialize empty dictionary to store data from each freq band
        self.multiFreqData = {}

        # loop over frequency band
        for fbSel,fbn in zip(fbandSel, fbands):
            all_data = []
            
            psel = []
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
                    if darkSel is None:
                        print "ERROR: no dark selection supplied"
                        return 0

                    fcmi, cmi, cmdti = getDarkModes(fdata, darkSel, [n_l,n_h],
                                                    df, nf, nsamps, par, tod)
                    fcm.append(fcmi)
                    cm.append(cmi)
                    cmdt.append(cmdti)

                r = lowFreqAnal(fdata, sel, [n_l,n_h], df, nsamps, scan_freq, par.get(parTag,{}), 
                                fcmodes = fcmi, respSel=respSel, flatfield=flatfield)
                
                psel.append(r["preSel"])
                corr.append(r["corr"])
                gain.append(np.abs(r["gain"]))
                norm.append(r["norm"])
                darkRatio.append(r["ratio"])
                
                if self._full:
                    all_data.append(r)
                    
            spsel = np.sum(psel,axis=0)
            Nmax = spsel.max()

            psel50 = spsel >= Nmax/2.
            
            for g,s in zip(gain,psel):
                g /= np.mean(g[psel50*s])
                
            gain = np.array(gain)
            gain[np.isnan(gain)] = 0.

            mgain = ma.MaskedArray(gain,~np.array(psel))
            mgain_mean = mgain.mean(axis=0)
            mcorr = ma.MaskedArray(corr,~np.array(psel))

            mcorr_max = mcorr.max(axis=0)
            mnorm = ma.MaskedArray(norm,~np.array(psel))
            mnorm_mean = mnorm.mean(axis=0)
                        
            if self._removeDark:
                mdarkRatio = ma.MaskedArray(darkRatio,~np.array(psel))
                mdarkRatio_mean = mdarkRatio.mean(axis=0)
                res['darkRatio'] = mdarkRatio_mean.data

            results = {
                "preLiveSel": psel50[fbSel],
                "liveSel": psel50[fbSel],
                "corrLive": mcorr_max.data[fbSel],
                "gainLive": mgain_mean.data[fbSel],
                "normLive": mnorm_mean.data[fbSel],
                "darkRatio": mdarkRatio_mean.data[fbSel]
            }

            if res.has_key('darkRatio'):
                self.crit["darkRatioLive"]["values"][fbSel] = res["darkRatio"][fbSel]

            multiFreqData[fbSel] = all_data
            
        # Undo flatfield correction
        self.crit["gainLive"] /= np.abs(ff)


    def lowFreqAnal(self, fdata, sel, frange, df, nsamps, scan_freq,
                    fcmodes=None, respSel=None, flatfield=None):
        """Find correlations and gains to the main common mode over a
        frequency range
        """
        # get relevant low freq data in the detectors selected
        lf_data = fdata[sel, frange[0]:frange[1]]
        ndet = len(sel)
        res = {}

        # Apply sine^2 taper to data
        if self._params.get("useTaper", False):
            taper = get_sine2_taper(frange, edge_factor = 6)
            lf_data *= np.repeat([taper],len(lf_data),axis=0)


        # Deproject correlated modes
        if fcmodes is not None:
            data_norm = np.linalg.norm(lf_data,axis=1)
            dark_coeff = []

            # actually do the deprojection here
            for m in fcmodes:
                coeff = numpy.dot(lf_data.conj(),m)
                lf_data -= numpy.outer(coeff.conj(),m)
                dark_coeff.append(coeff)

            # Reformat dark coefficients
            if len(dark_coeff) > 0:
                dcoeff = numpy.zeros([len(dark_coeff),ndet],dtype=complex)
                dcoeff[:,sel] = np.array(dark_coeff)

            # Get Ratio
            ratio = numpy.zeros(ndet,dtype=float)
            data_norm[data_norm==0.] = 1.
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
        ppar = self._params.get("presel",{})
        norm = np.zeros(ndet,dtype=float)
        fnorm = np.sqrt(np.abs(np.diag(c)))
        norm[sel] = fnorm*np.sqrt(2./nsamps)
        nnorm = norm/np.sqrt(nsamps)

        # get a range of valid norm values 
        nlim = ppar.get("normLimit",[0.,1e15])
        if np.ndim(nlim) == 0:
            nlim = [0, nlim]
        normSel = (nnorm > nlim[0])*(nnorm < nlim[1])
        
        # check which preselection is specified
        presel_method = ppar.get("method", "median")
        if presel_method is "median":
            sl = presel_by_median(cc, sel=normSel[sel], **presel_params)
            res["groups"] = None
            
        elif presel_method is "groups":
            G, ind, ld, smap = group_detectors(cc, sel=normSel[sel], **presel_params)
            sl = np.zeros(cc.shape[1], dtype=bool)
            sl[ld] = True
            res["groups"] = {
                "G": G,
                "ind": ind,
                "ld": ld,
                "smap": smap
            }
        else:
            raise "ERROR: Unknown preselection method"

        # The number of sels are just overwhelmingly confusing
        # To clarify for myself,
        # - normSel: selects the detectors with good norm
        # - sel: the initial selection of detectors specified
        #        for this case it is selection of dark detectors
        # - sl: is the preselected detectors from the median
        #       or group methods
        # Here it's trying to apply the preselection to the
        # dark selection
        preSel = sel.copy()
        preSel[sel] = sl


        # Apply gain ratio in case of multichroic
        if (flatfield is not None) and ("scale" in flatfield.fields):
            scl = flatfield.get_property("scale",det_uid=np.where(sel)[0],
                                         default = 1.)
            lf_data *= np.repeat([scl],lf_data.shape[1],axis=0).T

        # Get Correlations
        u, s, v = np.linalg.svd(lf_data[sl], full_matrices=False )
        corr = np.zeros(ndet)
        if par.get("doubleMode", False):
            corr[preSel] = np.sqrt(abs(u[:,0]*s[0])**2+abs(u[:,1]*s[1])**2)/fnorm[sl]
        else:
            corr[preSel] = np.abs(u[:,0])*s[0]/fnorm[sl]

        # Get Gains
        # data = CM * gain
        gain = np.zeros(ndet, dtype=complex)
        gain[preSel] = np.abs(u[:, 0])
        
        res.update({"preSel": preSel, "corr": corr, "gain": gain, "norm": norm, 
                    "dcoeff": dcoeff, "ratio": ratio, "cc": cc, "normSel": normSel})
        
        return res

    def getDarkModes(self, fdata, darkSel, frange, df, nf, nsamps, par, tod = None):
        """
        @brief Get dark or thermal modes from dark detectors and thermometer
               data.
        @return correlated modes in frequency and time domain, plus thermometer
                info.
        """
        n_l, n_h=frange
        fc_inputs = []

        # Dark detector drift
        if par["darkModesParams"].get("useDarks", False):
            dark_signal = fdata[darkSel,n_l:n_h].copy()
            fc_inputs.extend(list(dark_signal))

        # TEST CRYOSTAT TEMPERATURE
        if par["darkModesParams"].get("useTherm", False):
            thermometers = []
            for channel in par['thermParams']['channel']:
                # gather thermometer data
                thermometer = tod.get_hk( channel, fix_gaps=True)
                if len(np.diff(thermometer).nonzero()[0]) > 0:
                    thermometers.append(thermometer)
                    
            # perform a fourior analysis on thermometer data
            # and append to the dark detector data
            if len(thermometers) > 0:
                thermometers = np.array(thermometers)
                fth = np.fft.rfft( thermometers, nf )[:,n_l:n_h]
                fc_inputs.extend(list(fth))

        fc_inputs = np.array(fc_inputs)

        if par.get("useTaper",False):
            taper = get_sine2_taper(frange, edge_factor = 6)
            fc_inputs *= np.repeat([taper],len(fc_inputs),axis=0)

        # Normalize modes
        fc_inputs /= np.linalg.norm(fc_inputs, axis=1)[:, np.newaxis]

        # Obtain main svd modes to deproject from data
        if par["darkModesParams"].get("useSVD", False):
            Nmodes = par["darkModesParams"].get("Nmodes", None)
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


class GetSlowMode(Routine):
    def __init__(self, **params):
        self._params = params
        self._tod = params.get('tod', None)
        self._driftFilter = params.get('driftFilter', None)
        self._output_key = paramsg.et('output_key', None)

    def execute(self, store):
        tod = store.get(self._tod)

        # get the range of frequencies to work on
        n_l = 1
        n_h = nextregular(int(round(self._driftFilter/df))) + 1

        # extract the fourior modes
        lf_data = fdata[:,n_l:n_h]

        # calculate the mean fourior modes for live and dark 
        fcmL = lf_data[self.preLiveSel].mean(axis = 0)
        fcmD = lf_data[self.preDarkSel].mean(axis = 0)

        # get the common modes for both live and dark
        dsCM, dsCM_dt = get_time_domain_modes(fcmL,n_l, tod.nsamps, df)
        dsDCM, _ = get_time_domain_modes(fcmD,n_l, tod.nsamps, df)

        results = {
            "dsCM": dsCM,
            "dsDCM": dsDCM
        }
        store.set(self._output_key, results)

    
class Retrend(Routine):
    def __init__(self, **params):
        self._params = params


    def execute(self, store):
        # get the trend for the live detectors
        trL = numpy.array(trend).T[self.preLiveSel].mean(axis=0)
        trLt = trL[:,np.newaxis]

        # retrend the live detectors
        moby2.tod.retrend_tod(trLt, data = self.dsCM)
        
        # get the trend for the dark detectors
        trD = numpy.array(trend).T[self.preDarkSel].mean(axis=0)
        trDt = trD[:,np.newaxis]

        # retrend the dark detectors
        moby2.tod.retrend_tod(trDt, data = self.dsDCM)
        

class GetDriftErrors(Routine):
    def __init__(self, **params):
        self._params = params
        self._output_key = params.get('output_key', None)


    def execute(self, store):
        # Get Drift-Error
        DE = highFreqAnal(fdata, live, [n_l,n_h], self.ndata, nmodes = par["DEModes"], 
                          highOrder=False, preSel=self.preLiveSel)
        
        results = {
            "DELive": DE
        }

        store.set(self._output_key, results)


        
def highFreqAnal(fdata, sel, range, nsamps,  
                 nmodes=0, highOrder = False, preSel = None,
                 scanParams = None):
    """
    @brief Find noise RMS, skewness and kurtosis over a frequency band
    """
    ndet = len(sel)

    # get the high frequency fourior modes
    hf_data = fdata[sel,range[0]:range[1]]
    if nmodes > 0:
        # see if there is anything preselected
        if preSel is None:
            preSel = np.ones(sel.sum(),dtype=bool)
        else:
            preSel = preSel[sel]

        # find the correlation between different detectors
        c = np.dot(hf_data[preSel],hf_data[preSel].T.conjugate())

        # find the first few common modes in the detectors and
        # deproject them
        u, w, v = np.linalg.svd(c, full_matrices = 0)
        kernel = v[:nmodes]/np.repeat([np.sqrt(w[:nmodes])],len(c),axis=0).T
        modes = np.dot(kernel,hf_data[preSel])
        coeff = np.dot(modes,hf_data.T.conj())
        hf_data -= np.dot(coeff.T.conj(),modes)

    # compute the rms for the detectors
    rms = np.zeros(ndet)
    rms[sel] = np.sqrt(np.sum(abs(hf_data)**2,axis=1)/hf_data.shape[1]/nsamps)

    # if we are interested in high order effects, skew and kurtosis will
    # be calculated here
    if highOrder:
        hfd, _ = get_time_domain_modes( hf_data, 1, nsamps)
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


class AnalyzeMF(Routine):
    def __init__(self, **params):
        self._params = params
        self._midFreqFilter = params.get("midFreqFilter", None)
        self._output_key = params.get("output_key", None)


    def execute(self, store):
        # get the frequency range to work on
        n_l = int(self._midFreqFilter[0]/df))
        n_h = int(self._midFreqFilter[1]/df))        

        # perform a high frequency like analysis on this range
        MFE = highFreqAnal(fdata, live, [n_l,n_h], self.ndata, nmodes = par["MFEModes"], 
                             highOrder = False, preSel = self.preLiveSel)
        
        results = {
            "MFELive": MFE
        }

        store.set(self._output_key, results)


class AnalyzeHF(Routine):
    def __init__(self, **params):
        self._params = params

    def execute(self, store):
        # get the range of frequencies to work with
        n_l = int(round(par["highFreqFilter"][0]/df))
        n_h = int(round(par["highFreqFilter"][1]/df))
        # make sure that n_h is a number that's optimized in fft
        n_h = nextregular(n_h-n_l) + n_l

        # if partial is labeled the analysis will be carried out
        # for each individual scan between the turnning points
        if not(par["getPartial"]):
            rms, skewt, kurtt = highFreqAnal(fdata, live, [n_l,n_h], self.ndata, 
                                           nmodes=par["HFLiveModes"], 
                                           highOrder=True, preSel = self.preLiveSel)
        else:
            rms, skewt, kurtt, prms, pskewt, pkurtt=highFreqAnal(fdata, live, 
                                           [n_l,n_h], self.ndata, nmodes = par["HFLiveModes"], 
                                           highOrder=True, preSel=self.preLiveSel, 
                                           scanParams=self.scan)

            # store the statistics for partial
            results = {}
            results["partialRMSLive"] = np.zeros([self.ndet,self.chunkParams["N"]])
            results["partialSKEWLive"] = np.zeros([self.ndet,self.chunkParams["N"]]) 
            results["partialKURTLive"] = np.zeros([self.ndet,self.chunkParams["N"]]) 
            results["partialSKEWPLive"] = np.zeros([self.ndet,self.chunkParams["N"]])
            results["partialKURTPLive"] = np.zeros([self.ndet,self.chunkParams["N"]])
            results["partialRMSLive"][live] =  prms
            results["partialSKEWLive"][live] = pskewt.T[:,0]
            results["partialKURTLive"][live] = pkurtt.T[:,0]
            results["partialSKEWPLive"][live] = pskewt.T[:,1]
            results["partialKURTPLive"][live] = pkurtt.T[:,1]
            
        # store the statistics for global
        results["rmsLive"] = rms
        results["skewLive"] = np.zeros(self.ndet)
        results["kurtLive"] = np.zeros(self.ndet)
        results["skewpLive"] = np.zeros(self.ndet)
        results["kurtpLive"] = np.zeros(self.ndet)
        results["skewLive"][live] = skewt[0]
        results["kurtLive"][live] = kurtt[0]
        results["skewpLive"][live] = skewt[1]
        results["kurtpLive"][live] = kurtt[1]        

        # analyze the dark detectors for the same frequency range
        rms = highFreqAnal(fdata, dark, [n_l,n_h], self.ndata, nmodes = par["HFDarkModes"], 
                             highOrder = False, preSel = self.preDarkSel)
        results["rmsDark"] = rms


class AnalyzeAtm(Routine):
    def __init__(self, **params):
        """This routine does the atomosphere 1/f analysis"""
        self._params = params
        self._fitPowerLaw = params.get("fitPowerLaw", False)

    def execute(self, store):
        if self._fitPowerLaw:
            # look at both the live and dark detectors
            sel = preLiveSel + preDarkSel

            # combine the rms from both live and dark detectors
            rms = self.crit["rmsLive"]["values"] + self.crit["rmsDark"]["values"]

            # fit an atmosphere model with it 
            powLaw, level, knee = fit_atm(fdata, sel, dt, df, rms, self.scan_freq,
                                          **par.get("atmFit",{}))
            results = {
                "atmPowLaw": powLaw,
                "atmLevel": level,
                "atmKnee": knee
            }
            self.store(self._output_key, results)


            
def fit_atm(fdata, sel, dt, df, noise, scanf, 
            fminA=0.2, fmaxA=3., fmaxT = 10., 
            its = 1, width = 0.005):
    """
    Fit a power law to the atmosphere signal in a range of frequencies.
    """
    scale = 2*dt**2*df
    delta = 0.7
    kneeM = fmaxA+delta
    ind_ini = int(fminA/df)
    ps = np.power(np.abs(fdata[sel,:int(fmaxT/df)]),2)*scale
    # Get scan harmonics
    n_harm = int(np.ceil(fmaxT/scanf))
    i_harm = np.array(np.round(np.arange(1,n_harm+1)*scanf/df), dtype=int)
    i_harmw = i_harm.copy()
    di = int(width/df)
    for i in xrange(di):
        i_harmw = np.hstack([i_harmw,i_harm-(i+1)])
        i_harmw = np.hstack([i_harmw,i_harm+(i+1)])
    i_harmw.sort()
    # Iterate range fit
    for it in xrange(its):
        fmax = kneeM-delta-fminA
        imin = int(fminA/df); imax = int(fmax/df)
        psr = ps[:,imin:imax]
        log_ps = np.log(psr)
        freq = np.arange(imin,imax)*df
        log_freq = np.log(freq)
        w = np.diff(log_freq)
        w = np.hstack([w,w[-1]])
        iharm = i_harmw[(i_harmw>imin)*(i_harmw<imax)]-imin
        s = np.ones(len(freq),dtype=bool)
        s[iharm] = False
        m,n = np.polyfit(log_freq[s], log_ps[:,s].T, 1, w=w[s])
        pl = np.power(freq[np.newaxis,s].repeat(ps.shape[0],0),
                      m[:,np.newaxis].repeat(len(freq[s]),1))*\
                      np.exp(n[:,np.newaxis].repeat(len(freq[s]),1))
        c = np.sum(psr[:,s]*pl,axis=1)/np.sum(pl*pl,axis=1)
        level = np.exp(n)*c
        knee = np.power(noise[sel]/level,1./m)
        kneeM = np.median(knee)
    mA = np.zeros(fdata.shape[0])
    levelA = np.zeros(fdata.shape[0])
    kneeA = np.zeros(fdata.shape[0])
    mA[sel] = m
    levelA[sel] = level
    kneeA[sel] = knee
    return mA, levelA, kneeA            
