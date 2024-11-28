import numpy as np
import scipy.special as sp
import numexpr as ne

# there is a working example of the function call commented out
# in the bottom of this .py file

def biocal(Adet, sigdet, calcurve, yeartype, sar, bd, brok=0, abu=[], res=[]):
    """
    p95_4, p95, p68_2, calprob, medage = biocal(Adet, sigdet, calcurve, yeartype, sar, bd, brok=0, abu=[], res=[])

    Uses sediment (sediment accumulation rate, bioturbation) and abundance priors
    to calibrate a 14C date carried out on e.g. marine sediment cores.

    See also published paper: 
    B.C. Lougheed, 2022. Using sedimentological priors to improve 14C calibration 
    of bioturbated sediment archives. Radiocarbon, doi:10.1017/RDC.2021.116

    Parameters
    ----------

    Adet : float
        Laboratory radiocarbon determination in 14C years BP.
    
    sigdet : float
        Laboratory radiocarbon measurment error (1 sigma), in 14C years.
    
    calcurve : string
        Specify which calibration curve to choose (not case sensitive):

        'IntCal20', 'Marine20', 'SHCal20', 'IntCal13', 'Marine13',
        'SHCal13, 'IntCal09', 'Marine09', 'IntCal04', 'Marine04',
        'SHCal04, 'IntCal98', 'Marine98'

    yeartype : string
        Specifying how to handle years (not case sensitive):
        'Cal BP' or 'BCE/CE'
    
    sar : float
        Sediment accumulation rate (SAR) prior in cm/ka

    bd : float
        Bioturbation depth prior in cm.

    brok : float
        Optional. Estimated fragmented fraction of the sample's microfossil population that is fragmented,
        and therefore not picked. A value between 0 and 1. Will default to 0.

    abu : array-like
        Optional. Temporal abundance of the measured species. Two possibillities.

        First possibility (default):

        Empty, i.e. [], for temporally constant abundance.

        Second possiblity:

        Two coloumn array, i.e. [age abu] for temporally dynamic abundance.
        Each row contains an age (in years) in column 1 and relative
        abundance normalised to between 0 and 1 in column 2.
        Interpolation will be carried out to assign abundance to intermediate
        ages. Constant extrapolation of the edge values will be carried out to 
        assign abudance to ages outside the range provided.
        Ages (in years) should be entered either as Cal BP or BCE/CE,
        depending on what is set for yeartype.

    res : array-like
        Reservoir effect. Three possible types of input.

        First possibility (default):

        Empty, i.e. [], for no reservoir effect.

        Second possibility:

        Two column array, i.e. [reseff reserr], for a temporally constant
        reservoir effect (reseff) of reseff ± reserr. Reservoir effect
        to be entered as 14C yr offset from the desired calibration
        curve.
        
        Third possibility:

        Three column array, i.e. [age reseff reserr], for incorporating a temporally
        dynamic reservoir effect. Each row contains an age (in years) in 
        column 1 and reservoir effect in column 2 and the reservoir effect
        uncertainty in column 3. Reservoir effect should be entered in 14C yr 
        offset from the desired calibraiton curve and ages (in years) should 
        be entered either as Cal BP or BCE/CE, depending on what is set for yeartype.

    Returns
    -------

    p95_4 : ndarray
        Three column array containing 95.45% calibrated age probability
        range interval(s) calculated using highest posterior density.
        Each row contains a probability range in Cols 1 and 2, and
        the associated probability for that range in Col 3.
        Probabilities are normalised to between zero and one.

    p95 : ndarray
        Same as p95, but for 95% calibrated range.

    p68_2 : ndarray
        Same as p68_2, but for the 68.27% calibrated range.

    calprob : ndarray
        Two column array containing an annualised calibrated age
        probability density function for implementation in, e.g.,
        age modelling. Col 1 is a series of annual cal ages,
        Col 2 contains their associated probability. All probabilities
        are normalised such that they sum to 1.

    medage : float
        Median age calculated from calprob.


    Bryan Lougheed, Nov. 2024
    python 3.12.4, numpy 1.26.4
    """

    Adet = np.array(Adet)
    sigdet = np.array(sigdet)
    sar = np.array(sar)
    bd = np.array(bd)
    brok = np.array(brok)
    abu = np.array(abu)
    res = np.array(res)

    # Cal curve case and symbols
    headerlines = 11
    if calcurve.lower() == 'intcal20':
        calcurve = 'IntCal20'
        cite = '(Reimer et al., 2020)'
        curvetype = 'atm'
    elif calcurve.lower() == 'marine20':
        calcurve = 'Marine20'
        cite = '(Heaton et al., 2020)'
        curvetype = 'mar'
    elif calcurve.lower() == 'shcal20':
        calcurve = 'SHCal20'
        cite = '(Hogg et al., 2020)'
        curvetype = 'atm'
    elif calcurve.lower() == 'intcal13':
        calcurve = 'IntCal13'
        cite = '(Reimer et al., 2013)'
        curvetype = 'atm'
    elif calcurve.lower() == 'marine13':
        calcurve = 'Marine13'
        cite = '(Reimer et al., 2013)'
        curvetype = 'mar'
    elif calcurve.lower() == 'shcal13':
        calcurve = 'SHCal13'
        cite = '(Hogg et al., 2013)'
        curvetype = 'atm'
    elif calcurve.lower() == 'intcal09':
        calcurve = 'IntCal09'
        cite = '(Reimer et al., 2009)'
        curvetype = 'atm'
    elif calcurve.lower() == 'marine09':
        calcurve = 'Marine09'
        cite = '(Reimer et al., 2009)'
        curvetype = 'mar'
    elif calcurve.lower() == 'intcal04':
        calcurve = 'IntCal04'
        cite = '(Reimer et al., 2004)'
        curvetype = 'atm'
    elif calcurve.lower() == 'marine04':
        calcurve = 'Marine04'
        cite = '(Hughen et al., 2004)'
        curvetype = 'mar'
    elif calcurve.lower() == 'shcal04':
        calcurve = 'SHCal04'
        cite = '(McCormac et al., 2004)'
        curvetype = 'atm'
    elif calcurve.lower() == 'intcal98':
        headerlines = 18
        calcurve = 'IntCal98'
        cite = '(Stuiver et al., 1998)'
        curvetype = 'atm'
    elif calcurve.lower() == 'marine98':
        headerlines = 18
        calcurve = 'Marine98'
        cite = '(Stuiver et al., 1998)'
        curvetype = 'mar'
    else:
        raise ValueError('Calibration curve '+calcurve+' unknown. Please specify a valid calibration curve (see help for options)')

    # Process some of the user inputs
    sar = sar / 1000  # sar from cm/ka to cm/a

    # Load cal curve data into workspace
    try:
        contents = np.loadtxt('private/'+calcurve+'.14c', delimiter=',', skiprows=headerlines)
    except OSError:
        raise FileNotFoundError('Could not find '+calcurve+'.14c. Typo?')
    curvecal = np.flipud(contents[:, 0])
    curve14c = np.flipud(contents[:, 1])
    curveerr = np.flipud(contents[:, 2])

    # Process reservoir effect
    if res.size == 0:  # no reservoir effect
        pass # do nothing
    elif res.size == 2:  # constant reservoir effect. Still in 14C yr, so add to cal curve
        curve14c = curve14c + res[0]
        curveerr = np.sqrt(curveerr**2 + res[1]**2)
    elif res.shape[1] == 3:  # temporally dynamic reservoir effect
        try:
            if yeartype.lower() == 'bce/ce':
                res[:, 0] = (1950 - res[:, 0])
            res = np.sort(res, axis=0)
            reseff = np.interp(curvecal, res[:, 0], res[:, 1])
            reseff[curvecal < res[0, 0]] = res[0, 1]  # extrapolate using constant value
            reseff[curvecal > res[-1, 0]] = res[-1, 1]  # extrapolate using constant value
            curve14c = curve14c + reseff
            reserr = np.interp(curvecal, res[:, 0], res[:, 2])
            reserr[curvecal < res[0, 0]] = res[0, 2]  # extrapolate using constant value
            reserr[curvecal > res[-1, 0]] = res[-1, 2]  # extrapolate using constant value
            curveerr = np.sqrt(curveerr ** 2 + reserr ** 2)
        except:
            raise ValueError('Reservoir effect not entered correctly')
    else:
        raise ValueError('Reservoir effect not entered correctly')

    # linearly interpolate cal curve to 1 cal year resolution
    calres = 1
    curve14c = np.interp(np.arange(min(curvecal), max(curvecal) + 1, calres), curvecal, curve14c)
    curveerr = np.interp(np.arange(min(curvecal), max(curvecal) + 1, calres), curvecal, curveerr)
    curvecal = np.arange(min(curvecal), max(curvecal) + 1, calres)
    # convert interpolated curve to F14C activity
    curvef14 = np.exp(curve14c / -8033)
    curveerr = curvef14 * curveerr / 8033

    # collect abundance information
    if abu.size == 0:
        abus = 1
    elif abu.shape[1] == 2:
        try:
            if yeartype.lower() == 'bce/ce':
                abu[:,0] = (1950 - abu[:,0])
            abus = np.interp(curvecal, abu[:,0], abu[:,1])
            abus[curvecal < abu[0,0]] = abu[0,1]  # extrapolate using constant value
            abus[curvecal > abu[-1,0]] = abu[-1,1]  # extrapolate using constant value
        except:
            raise ValueError('Abundance not entered correctly')
    else:
        raise ValueError('Abundance not entered correctly')

    # make prior distribution (ppri) of relative age for the sample
    rnga = np.arange(0, round((bd * 5 / sar)) + 1)  # discrete relative age range of prior (5 bioturbation depths)
    ppri = np.exp(-(rnga * sar) / bd)  # Eq. 1 in manuscript
    if brok > 0:
        rk = (-bd * np.log(brok)) / sar  # Eq. 3 in manuscript
        ppri = ppri[rnga < rk]  # trim prior distribution to whole microfossils only
    ppri /= np.sum(ppri)  # normalise prior dist such that it sums to 1

    # start (ts) and end (te) of sliding windows, based on C14 errors and ppri length
    # get estimate (in the future, when computers get faster and have more memory, we can just do
    # the entire curve and not need to trim here)
    ts = np.argmin(np.abs(curve14c - (Adet - (3*sigdet) - ppri.size*0.75)))
    if ts.size > 1:
        ts = ts[0]
    te = np.argmin(np.abs(curve14c - (Adet + (3*sigdet) + ppri.size*0.75)))
    if te.size > 1:
        te = te[-1]

    # now also convert Adet and sigdet to F14C activity
    Adet = np.exp(Adet / -8033)
    sigdet = Adet * sigdet / 8033

    # initiate prob matrix pmat for all windows t, with space for ppri tail at final t
    pmat = np.zeros((te - ts + 1, te-ts+ppri.size))  #

    # trim cal curve to match prob matrix
    ind = (curvecal >= curvecal[ts]) & (curvecal <= curvecal[ts] + pmat.shape[1] - 1)
    curvef14 = curvef14[ind]
    curveerr = curveerr[ind]
    curvecal = curvecal[ind]
    # also trim abus
    if abu.size > 0:
        abus = abus[ind]

    # retrim pmat to calcurve (i.e. in case it exceeds end of cal curve)
    pmat = pmat[:, 0:curvecal.size]

    # populate pmat with ppri. each row (t) of pmat contains ppri placed at a new sliding window starting at t
    # where ppri(1) is placed at t
    S = pmat.shape[1]
    N = ppri.size

    for t in range(pmat.shape[0]):  # could perhaps be vectorised, but very fast as it is
        if S-t+1 >= N:  # enough space for the entire ppri
            pmat[t, t:t+N] = ppri
        else:  # nearing end of calcurve (and thus pmat)
            pmat[t, t:] = ppri[0 : S-t]
    pmat /= np.sum(pmat, axis=1, keepdims=True)  # normalise all rows of pmat
    pmat *= abus  # multiply by abundance probability
    pmat /= np.sum(pmat, axis=1, keepdims=True)  # normalise again

    # implicitly prep some stuff for vectorisation in numpy
    curvef14 = curvef14.reshape(1, curvef14.size)
    curveerr = curveerr.reshape(1, curveerr.size)

    ### Equation 6 in manuscript: calculate p14c(T|t) 
    # Fully vectorised here, each row corresponds to each instance of t, the we sum the rows
    # Line below fully vectorised, it works, but still slow and memory hogging in numpy:
    # p14cTt = np.sum(1 / (curveerr * np.sqrt(2*np.pi)) * np.exp(-(curvef14 - curvef14.T) ** 2 / (2*curveerr**2)), axis=1)
    # Below a faster/uglier way using numexpr, also less memory hog. Almost as good as Matlab
    diffmat = curvef14 - curvef14.T
    inv_curveerr = 1 / (curveerr * np.sqrt(2*np.pi))
    p14cTt = ne.evaluate('inv_curveerr * exp(-diffmat**2 / (2*curveerr**2))')
    p14cTt = np.sum(p14cTt, axis=1)

    ### Equation 7, hdet(t). Once again vectorised, each row of hdet corresponds to each instance of t.
    p14cTt = pmat * p14cTt[np.newaxis,:] # First part of Eq 7. (except abundance already taken care of above in line 263)
    p14cTt /= np.sum(p14cTt, axis=1, keepdims=True)  # normalise all rows. Each row is each instance t. (keep dims for vectorisation)
    hdet = np.sum(curvef14*p14cTt, axis=1)  # the final part of Eq. 7

    ### Equation 8 probability for each sliding window t based on the closeness of its hdet(t)
    # to labdet as calculated using normal pdf of Adet and sigdet. Vectorised.
    phdet = 1 / (sigdet * np.sqrt(2*np.pi)) * np.exp(-(hdet - Adet) ** 2 / (2 * sigdet**2))

    # Construct final cal age probability distribution
    pmat *= phdet[:, np.newaxis]  # Mcal(T)
    pcal = np.sum(pmat, axis=0)  # Pcal(T)
    pcal /= np.sum(pcal)  # normalise

    # output calprob to user
    calprob = np.column_stack((curvecal, pcal))

    # calculate median cal age
    median_ind = np.argmin(np.abs(np.cumsum(calprob[:,1]) - 0.5))
    medage = np.round(np.mean(calprob[median_ind, 0]))  # in case more than one

    # find 68.2% and 95.4% cal age credible intervals using highest posterior density (HPD)
    def hpdcalc(values, probs, valres, siglevel):
        """
        Calculates highest posterior density (HPD) credible interval.

        Parameters
        ----------
        values : array-like
            Vector of values, equally spaced, monotonic
        probs : array-like
            Vector of probabilities for each value. Total probability should sum to 1.
        valres : float
            The spacing in values, e.g., 1 for 1 year (if the values are years)
        siglevel : float
            Desired significance level, e.g., 0.95 for 95%

        Returns
        -------
        pranges : ndarray
            n by 3 matrix. Value interval(s) calculated using highest posterior density.
            Each row contains a probability range in Cols 1 and 2, and
            the associated probability for that range in Col 3.
            Probabilities are normalized to between zero and one.

        B.C. Lougheed
        """
        if not (isinstance(values, np.ndarray) and isinstance(probs, np.ndarray)):
            raise ValueError('Check that values and probs are both vectors')
        if len(values) != len(probs):
            raise ValueError('Check that values and probs have same number of elements')

        hpd = np.column_stack( (values, probs, np.full_like(values,np.nan)) )
        hpd = hpd[np.argsort(hpd[:,1])]
        hpd[:,2] = np.cumsum(hpd[:,1])
        hpd = hpd[hpd[:,2] >= 1 - siglevel]
        hpd = hpd[np.argsort(hpd[:,0])]
        ind1 = np.where(np.diff(hpd[:,0]) > valres + valres * 0.01)[0]  # 0.01 to deal with rounding errors

        if ind1.size == 0:
            pranges = np.array([[hpd[-1, 0], hpd[0, 0], np.sum(hpd[:, 1])]])
        else:
            ind2 = np.concatenate(([0], ind1, ind1 + 1, [len(hpd) - 1]))
            pranges = np.empty((len(ind1) + 1, 3))
            # need to double check this
            for i in range(len(ind1) + 1):
                pranges[i, 0] = hpd[ind2[2 * i + 1], 0]
                pranges[i, 1] = hpd[ind2[2 * i], 0]
                pranges[i, 2] = np.sum(hpd[ind2[2 * i]:ind2[2 * i + 1] + 1, 1])
            pranges = np.flipud(pranges)
        return pranges
    
    p95_4 = hpdcalc(calprob[:,0], calprob[:,1], 1, sp.erf(2 / np.sqrt(2)))
    p95 = hpdcalc(calprob[:,0], calprob[:,1], 1, 0.95)
    p68_2 = hpdcalc(calprob[:,0], calprob[:,1], 1, sp.erf(1 / np.sqrt(2)))

    # convert output to BCE/CE if requested
    if yeartype.lower() == 'bce/ce':
        medage = (medage - 1950) * -1
        calprob[:,0] = (calprob[:,0] - 1950) * -1
        p95_4[:,0:1] = (p95_4[:,0:1] - 1950) * -1
        p95[:,0:1] = (p95[:,0:1] - 1950) * -1
        p68_2[:,0:1] = (p68_2[:,0:1] - 1950) * -1

    return p95_4, p95, p68_2, calprob, medage

# # Working Example
# Adet = 12500
# sigdet= 40
# calcurve = 'Marine20'
# yeartype = 'Cal BP'
# sar = 5
# bd = 8
# brok = 0.1
# abu = []
# res = [650, 40]

# import time as time
# tstart = time.time()
# p95_4, p95, p68_2, calprob, medage = biocal(Adet, sigdet, calcurve, yeartype, sar, bd, brok, abu, res)
# print(str(round(time.time()-tstart,2))+' seconds')

# import matplotlib.pyplot as plt
# plt.plot(calprob[:,0],calprob[:,1])
# plt.xlabel('Cal yr BP')
# plt.ylabel('Probability')
# plt.show()