from __future__ import division
import numpy as np

def nextregular(n):
    while not checksize(n): n+=1
    return n

def checksize(n):
    while not (n%16): n//=16
    while not (n%13): n//=13
    while not (n%11): n//=11
    while not (n%9): n//=9
    while not (n%7): n//=7
    while not (n%5): n//=5
    while not (n%3): n//=3
    while not (n%2): n//=2
    return (1 if n == 1 else 0)


def presel_by_median(cc, sel=None, **kwargs):
    """
    minCorr: minimum correlation requiered for preselection
    superMinCorr: minimum correlation requiered in case minCorr produces less than
        max(<<minSel>>,numberOfDetectors/<<minFrac>> detectors
    minSel: minimum number of detectors preselected
    minFrac: determines the minimum number of detectors preselected by determining a 
        fraction of the number of detectors available.
    Note: to go back to c9 you can set:
        superMinCorr = 0.5
        minSel = 0
        minFrac = 10000
    """
    if sel is None:
        sel = np.ones(cc.shape[0],dtype=bool)
        
    minCorr = kwargs.get("minCorr", 0.6)
    superMinCorr = kwargs.get("superMinCorr", 0.3)
    minSel = kwargs.get("minSel", 10)
    minFrac = kwargs.get("minFrac", 10)

    # select those detectors whose medium are above a specified threshold
    sl = (np.median(abs(cc),axis=0) > minCorr)*sel
    
    if kwargs.get("forceSel") is not None:
        sl *= kwargs.get("forceSel") # NOT PRETTY
        
    if sl.sum() < np.max([cc.shape[0]/minFrac,minSel]):
        print("ERROR: did not find any valid detectors for low frequency analysis.")
        sl = (np.median(abs(cc),axis=0) > superMinCorr)*sel
        if sl.sum() < minSel:
            raise RuntimeError("PRESELECTION FAILED, did not find any valid detectors for low frequency analysis.")
    else:
        sl = ((abs(cc[sl]).mean(axis=0)-1./len(cc[sl]))*len(cc[sl])/(len(cc[sl])-1) > minCorr)*sel
    return sl


def group_detectors(cc, sel = None, **kwargs):
    """
    Groups detectors according to their correlation.
    Returns:
        G: list of lists of detectors containing the indexes of the detectors in each group
        ind: index of the last group included in the live detector preselection
        ld: indexes of detectors from the main correlated groups
    Note: Indexes are provided according to the correlation matrix given
    """
    thr0 = kwargs.get("initCorr",0.99)
    thr1 = kwargs.get("groupCorr",0.8)
    thrg = kwargs.get("minCorr",0.6)
    dthr = kwargs.get("deltaCorr", 0.005)
    Nmin = kwargs.get("Nmin",20)
    Gmax = kwargs.get("Gmax",5)

    if sel is None: sel = np.ones(cc.shape[0],dtype=bool)
    smap = np.where(sel)[0]
    scc = cc[sel][:,sel]

    G = []
    g0 = []
    allind = np.arange(scc.shape[0])
    ss = np.zeros(scc.shape[0],dtype=bool)
    thr = thr0
    while ss.sum() < len(allind):
        if np.sum(~ss) <= Nmin or len(G) >= Gmax: 
            G.append(smap[np.where(~ss)[0]])
            break
 
        ind = allind[~ss]
        N = np.sum(~ss)
        cco = scc[~ss][:,~ss]
 
        # Find reference mode
        n0 = np.sum(np.abs(cco)>thr,axis=0)
        imax = np.argmax(n0)
        if n0[imax] < np.min([Nmin,N/2]):
            thr -= dthr
            continue
 
        # Find initial set of strongly correlated modes
        gg = np.where(np.abs(cco[imax])>thr)[0]
        s = np.argsort(cco[imax][gg])
        g = ind[gg[s]].tolist()
 
        # Extend set until thr1
        while thr > thr1:
            thr -= dthr
            sg = np.ones(scc.shape[0],dtype=bool)
            sg[g] = False
            sg[ss] = False
            ind = np.where(sg)[0]
  
            if np.sum(sg) <= Nmin:
                g0.extend(np.where(sg)[0])
                break
  
            cci = scc[g][:,sg]
            g0 = np.where(~np.any(np.abs(cci)<thr,axis=0))[0].tolist()
            g.extend(ind[g0])
 
        # Append new group result
        G.append(smap[g])
        ss[g] = True
        #print len(g), thr
        thr = thr0

    ind = 0
    ld = G[ind].tolist()
    while (ind < len(G)-1):
        if np.mean(np.abs(cc[G[0]][:,G[ind+1]])) > thrg:
            ind += 1
            ld.extend(G[ind].tolist())
        else:
            break

    return G, ind, ld, smap    


def get_sine2_taper(frange, edge_factor = 6):
    # Generate a frequency space taper to reduce ringing in lowFreqAnal
    band = frange[1]-frange[0]
    edge = band/edge_factor
    x = np.arange(edge, dtype=float) / (edge-1)
    taper = np.ones(band)
    taper[:edge] = np.sin(x*np.pi/2)**2
    taper[-edge:] = np.sin((x+1)*np.pi/2)**2
    return taper

def get_iharm(frange, df, scan_freq, wide = False):
    # Get the harmonic mode of the scan frequency
    n_harm = int(np.ceil(frange[1]*df/scan_freq))
    f_harm = (np.arange(n_harm)+1)*scan_freq
    if wide:
        i_harm = np.array(np.sort(np.hstack(
                      [np.round(f_harm/df-1),
                       np.round(f_harm/df),
                       np.round(f_harm/df+1)])),dtype=int)
    else:
        i_harm = np.array(np.round(f_harm/df), dtype=int)
    i_harm = i_harm[(i_harm>=frange[0])*(i_harm<frange[1])] - frange[0]
    return i_harm


def get_time_domain_modes(fmodes, n_l, nsamps, df=1.,):
    # Get modes in time domain
    if fmodes.ndim == 1:
        fmodes = fmodes[np.newaxis,:]
    fcm = np.hstack([np.zeros((len(fmodes),n_l)),
                        fmodes[:,:-1],
                        np.expand_dims(np.real(fmodes[:,-1]),1)]) 
    modes = np.fft.irfft(fcm)
    modes_dt = 1./modes.shape[1]/df
    modes *= np.sqrt(2.*fmodes.shape[1]/nsamps)
    return modes, modes_dt
