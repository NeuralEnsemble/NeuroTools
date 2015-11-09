"""
NeuroTools.analysis
===================

A collection of analysis functions that may be used by NeuroTools.signals or other packages.

.. currentmodule:: NeuroTools.analysis

Classes
-------
.. autosummary::

   TuningCurve

Functions
---------
.. autosummary::
   :nosignatures:

   ccf
   crosscorrelate
   make_kernel
   simple_frequency_spectrum

"""

import numpy as np

from NeuroTools import check_dependency

HAVE_MATPLOTLIB = check_dependency('matplotlib')
if HAVE_MATPLOTLIB:
    import matplotlib
    matplotlib.use('Agg')
else:
    MATPLOTLIB_ERROR = "The matplotlib package was not detected"

HAVE_PYLAB = check_dependency('pylab')
if HAVE_PYLAB:
    import pylab
else:
    PYLAB_ERROR = "The pylab package was not detected"


def ccf(x, y, axis=None):
    """Fast cross correlation function based on fft.

    Computes the cross-correlation function of two series.
    Note that the computations are performed on anomalies (deviations from
    average).
    Returns the values of the cross-correlation at different lags.

    Parameters
    ----------
    x, y : 1D MaskedArrays
        The two input arrays.
    axis : integer, optional
        Axis along which to compute (0 for rows, 1 for cols).
        If `None`, the array is flattened first.

    Examples
    --------
    >>> z = arange(5)
    >>> ccf(z,z)
    array([  3.90798505e-16,  -4.00000000e-01,  -4.00000000e-01,
            -1.00000000e-01,   4.00000000e-01,   1.00000000e+00,
             4.00000000e-01,  -1.00000000e-01,  -4.00000000e-01,
            -4.00000000e-01])

    """
    assert x.ndim == y.ndim, "Inconsistent shape !"
#    assert(x.shape == y.shape, "Inconsistent shape !")
    if axis is None:
        if x.ndim > 1:
            x = x.ravel()
            y = y.ravel()
        npad = x.size + y.size
        xanom = (x - x.mean(axis=None))
        yanom = (y - y.mean(axis=None))
        Fx = np.fft.fft(xanom, npad, )
        Fy = np.fft.fft(yanom, npad, )
        iFxy = np.fft.ifft(Fx.conj() * Fy).real
        varxy = np.sqrt(np.inner(xanom, xanom) * np.inner(yanom, yanom))
    else:
        npad = x.shape[axis] + y.shape[axis]
        if axis == 1:
            if x.shape[0] != y.shape[0]:
                raise ValueError("Arrays should have the same length!")
            xanom = (x - x.mean(axis=1)[:, None])
            yanom = (y - y.mean(axis=1)[:, None])
            varxy = np.sqrt((xanom * xanom).sum(1) *
                            (yanom * yanom).sum(1))[:, None]
        else:
            if x.shape[1] != y.shape[1]:
                raise ValueError("Arrays should have the same width!")
            xanom = (x - x.mean(axis=0))
            yanom = (y - y.mean(axis=0))
            varxy = np.sqrt((xanom * xanom).sum(0) * (yanom * yanom).sum(0))
        Fx = np.fft.fft(xanom, npad, axis=axis)
        Fy = np.fft.fft(yanom, npad, axis=axis)
        iFxy = np.fft.ifft(Fx.conj() * Fy, n=npad, axis=axis).real
    # We just turn the lags into correct positions:
    iFxy = np.concatenate((iFxy[len(iFxy) / 2:len(iFxy)],
                           iFxy[0:len(iFxy) / 2]))
    return iFxy / varxy

from NeuroTools.plotting import get_display, set_labels

HAVE_PYLAB = check_dependency('pylab')


def crosscorrelate(sua1, sua2, lag=None, n_pred=1, predictor=None,
                   display=False, kwargs={}):
    """Cross-correlation between two series of discrete events (e.g. spikes).

    Calculates the cross-correlation between
    two vectors containing event times.
    Returns ``(differeces, pred, norm)``. See below for details.

    Adapted from original script written by Martin P. Nawrot for the
    FIND MATLAB toolbox [1]_.

    Parameters
    ----------
    sua1, sua2 : 1D row or column `ndarray` or `SpikeTrain`
        Event times. If sua2 == sua1, the result is the autocorrelogram.
    lag : float
        Lag for which relative event timing is considered
        with a max difference of +/- lag. A default lag is computed
        from the inter-event interval of the longer of the two sua
        arrays.
    n_pred : int
        Number of surrogate compilations for the predictor. This
        influences the total length of the predictor output array
    predictor : {None, 'shuffle'}
        Determines the type of bootstrap predictor to be used.
        'shuffle' shuffles interevent intervals of the longer input array
        and calculates relative differences with the shorter input array.
        `n_pred` determines the number of repeated shufflings, resulting
        differences are pooled from all repeated shufflings.
    display : boolean
        If True the corresponding plots will be displayed. If False,
        int, int_ and norm will be returned.
    kwargs : dict
        Arguments to be passed to np.histogram.

    Returns
    -------
    differences : np array
        Accumulated differences of events in `sua1` minus the events in
        `sua2`. Thus positive values relate to events of `sua2` that
        lead events of `sua1`. Units are the same as the input arrays.
    pred : np array
        Accumulated differences based on the prediction method.
        The length of `pred` is ``n_pred * length(differences)``. Units are
        the same as the input arrays.
    norm : float
        Normalization factor used to scale the bin heights in `differences` and
        `pred`. ``differences/norm`` and ``pred/norm`` correspond to the linear
        correlation coefficient.

    Examples
    --------
    >> crosscorrelate(np_array1, np_array2)
    >> crosscorrelate(spike_train1, spike_train2)
    >> crosscorrelate(spike_train1, spike_train2, lag = 150.0)
    >> crosscorrelate(spike_train1, spike_train2, display=True,
                      kwargs={'bins':100})

    See also
    --------
    ccf

    .. [1] Meier R, Egert U, Aertsen A, Nawrot MP, "FIND - a unified framework
       for neural data analysis"; Neural Netw. 2008 Oct; 21(8):1085-93.

    """
    assert predictor is 'shuffle' or predictor is None, "predictor must be \
    either None or 'shuffle'. Other predictors are not yet implemented."

    #Check whether sua1 and sua2 are SpikeTrains or arrays
    sua = []
    for x in (sua1, sua2):
        #if isinstance(x, SpikeTrain):
        if hasattr(x, 'spike_times'):
            sua.append(x.spike_times)
        elif x.ndim == 1:
            sua.append(x)
        elif x.ndim == 2 and (x.shape[0] == 1 or x.shape[1] == 1):
            sua.append(x.ravel())
        else:
            raise TypeError("sua1 and sua2 must be either instances of the" \
                            "SpikeTrain class or column/row vectors")
    sua1 = sua[0]
    sua2 = sua[1]

    if sua1.size < sua2.size:
        if lag is None:
            lag = np.ceil(10*np.mean(np.diff(sua1)))
        reverse = False
    else:
        if lag is None:
            lag = np.ceil(20*np.mean(np.diff(sua2)))
        sua1, sua2 = sua2, sua1
        reverse = True

    #construct predictor
    if predictor is 'shuffle':
        isi = np.diff(sua2)
        sua2_ = np.array([])
        for ni in xrange(1,n_pred+1):
            idx = np.random.permutation(isi.size-1)
            sua2_ = np.append(sua2_, np.add(np.insert(
                (np.cumsum(isi[idx])), 0, 0), sua2.min() + (
                np.random.exponential(isi.mean()))))

    #calculate cross differences in spike times
    differences = np.array([])
    pred = np.array([])
    for k in xrange(0, sua1.size):
        differences = np.append(differences, sua1[k] - sua2[np.nonzero(
            (sua2 > sua1[k] - lag) & (sua2 < sua1[k] + lag))])
    if predictor == 'shuffle':
        for k in xrange(0, sua1.size):
            pred = np.append(pred, sua1[k] - sua2_[np.nonzero(
                (sua2_ > sua1[k] - lag) & (sua2_ < sua1[k] + lag))])
    if reverse is True:
        differences = -differences
        pred = -pred

    norm = np.sqrt(sua1.size * sua2.size)

    # Plot the results if display=True
    if display:
        subplot = get_display(display)
        if not subplot or not HAVE_PYLAB:
            return differences, pred, norm
        else:
            # Plot the cross-correlation
            try:
                counts, bin_edges = np.histogram(differences, **kwargs)
                edge_distances = np.diff(bin_edges)
                bin_centers = bin_edges[1:] - edge_distances/2
                counts = counts / norm
                xlabel = "Time"
                ylabel = "Cross-correlation coefficient"
                #NOTE: the x axis corresponds to the upper edge of each bin
                subplot.plot(bin_centers, counts, label='cross-correlation', color='b')
                if predictor is None:
                    set_labels(subplot, xlabel, ylabel)
                    pylab.draw()
                elif predictor is 'shuffle':
                    # Plot the predictor
                    norm_ = norm * n_pred
                    counts_, bin_edges_ = np.histogram(pred, **kwargs)
                    counts_ = counts_ / norm_
                    subplot.plot(bin_edges_[1:], counts_, label='predictor')
                    subplot.legend()
                    pylab.draw()
            except ValueError:
                print("There are no correlated events within the selected lag"\
                " window of %s" % lag)
    else:
        return differences, pred, norm


def _dict_max(D):
    """For a dict containing numerical values, return the key for the
    highest value. If there is more than one item with the same highest
    value, return one of them (arbitrary - depends on the order produced
    by the iterator).

    """
    max_val = max(D.values())
    for k in D:
        if D[k] == max_val:
            return k


def make_kernel(form, sigma, time_stamp_resolution, direction=1):
    """Creates kernel functions for convolution.

    Constructs a numeric linear convolution kernel of basic shape to be used
    for data smoothing (linear low pass filtering) and firing rate estimation
    from single trial or trial-averaged spike trains.

    Exponential and alpha kernels may also be used to represent postynaptic
    currents / potentials in a linear (current-based) model.

    Adapted from original script written by Martin P. Nawrot for the
    FIND MATLAB toolbox [1]_ [2]_.

    Parameters
    ----------
    form : {'BOX', 'TRI', 'GAU', 'EPA', 'EXP', 'ALP'}
        Kernel form. Currently implemented forms are BOX (boxcar),
        TRI (triangle), GAU (gaussian), EPA (epanechnikov), EXP (exponential),
        ALP (alpha function). EXP and ALP are aymmetric kernel forms and
        assume optional parameter `direction`.
    sigma : float
        Standard deviation of the distribution associated with kernel shape.
        This parameter defines the time resolution (in ms) of the kernel estimate
        and makes different kernels comparable (cf. [1] for symetric kernels).
        This is used here as an alternative definition to the cut-off
        frequency of the associated linear filter.
    time_stamp_resolution : float
        Temporal resolution of input and output in ms.
    direction : {-1, 1}
        Asymmetric kernels have two possible directions.
        The values are -1 or 1, default is 1. The
        definition here is that for direction = 1 the
        kernel represents the impulse response function
        of the linear filter. Default value is 1.

    Returns
    -------
    kernel : array_like
        Array of kernel. The length of this array is always an odd
        number to represent symmetric kernels such that the center bin
        coincides with the median of the numeric array, i.e for a
        triangle, the maximum will be at the center bin with equal
        number of bins to the right and to the left.
   norm : float
        For rate estimates. The kernel vector is normalized such that
        the sum of all entries equals unity sum(kernel)=1. When
        estimating rate functions from discrete spike data (0/1) the
        additional parameter `norm` allows for the normalization to
        rate in spikes per second.

        For example:
        ``rate = norm * scipy.signal.lfilter(kernel, 1, spike_data)``
    m_idx : int
        Index of the numerically determined median (center of gravity)
        of the kernel function.

    Examples
    --------
    To obtain single trial rate function of trial one should use::

        r = norm * scipy.signal.fftconvolve(sua, kernel)

    To obtain trial-averaged spike train one should use::

        r_avg = norm * scipy.signal.fftconvolve(sua, np.mean(X,1))

    where `X` is an array of shape `(l,n)`, `n` is the number of trials and
    `l` is the length of each trial.

    See also
    --------
    SpikeTrain.instantaneous_rate
    SpikeList.averaged_instantaneous_rate

    .. [1] Meier R, Egert U, Aertsen A, Nawrot MP, "FIND - a unified framework
       for neural data analysis"; Neural Netw. 2008 Oct; 21(8):1085-93.

    .. [2] Nawrot M, Aertsen A, Rotter S, "Single-trial estimation of neuronal
       firing rates - from single neuron spike trains to population activity";
       J. Neurosci Meth 94: 81-92; 1999.

    """
    assert form.upper() in ('BOX','TRI','GAU','EPA','EXP','ALP'), "form must \
    be one of either 'BOX','TRI','GAU','EPA','EXP' or 'ALP'!"

    assert direction in (1,-1), "direction must be either 1 or -1"

    SI_sigma = sigma / 1000. #convert to SI units (ms -> s)

    SI_time_stamp_resolution = time_stamp_resolution / 1000. #convert to SI units (ms -> s)

    norm = 1./SI_time_stamp_resolution

    if form.upper() == 'BOX':
        w = 2.0 * SI_sigma * np.sqrt(3)
        width = 2 * np.floor(w / 2.0 / SI_time_stamp_resolution) + 1  # always odd number of bins
        height = 1. / width
        kernel = np.ones((1, width)) * height  # area = 1

    elif form.upper() == 'TRI':
        w = 2 * SI_sigma * np.sqrt(6)
        halfwidth = np.floor(w / 2.0 / SI_time_stamp_resolution)
        trileft = np.arange(1, halfwidth + 2)
        triright = np.arange(halfwidth, 0, -1)  # odd number of bins
        triangle = np.append(trileft, triright)
        kernel = triangle / triangle.sum()  # area = 1

    elif form.upper() == 'EPA':
        w = 2.0 * SI_sigma * np.sqrt(5)
        halfwidth = np.floor(w / 2.0 / SI_time_stamp_resolution)
        base = np.arange(-halfwidth, halfwidth + 1)
        parabula = base**2
        epanech = parabula.max() - parabula  # inverse parabula
        kernel = epanech / epanech.sum()  # area = 1

    elif form.upper() == 'GAU':
        w = 2.0 * SI_sigma * 2.7  # > 99% of distribution weight
        halfwidth = np.floor(w / 2.0 / SI_time_stamp_resolution)  # always odd
        base = np.arange(-halfwidth, halfwidth + 1) * SI_time_stamp_resolution
        g = np.exp(-(base**2) / 2.0 / SI_sigma**2) / SI_sigma / np.sqrt(2.0 * np.pi)
        kernel = g / g.sum()

    elif form.upper() == 'ALP':
        w = 5.0 * SI_sigma
        alpha = np.arange(1, (2.0 * np.floor(w / SI_time_stamp_resolution / 2.0) + 1) + 1) * SI_time_stamp_resolution
        alpha = (2.0 / SI_sigma**2) * alpha * np.exp(-alpha * np.sqrt(2) / SI_sigma)
        kernel = alpha / alpha.sum()  # normalization
        if direction == -1:
            kernel = np.flipud(kernel)

    elif form.upper() == 'EXP':
        w = 5.0 * SI_sigma
        expo = np.arange(1, (2.0 * np.floor(w / SI_time_stamp_resolution / 2.0) + 1) + 1) * SI_time_stamp_resolution
        expo = np.exp(-expo / SI_sigma)
        kernel = expo / expo.sum()
        if direction == -1:
            kernel = np.flipud(kernel)

    kernel = kernel.ravel()
    m_idx = np.nonzero(kernel.cumsum() >= 0.5)[0].min()

    return kernel, norm, m_idx


def simple_frequency_spectrum(x):
    """Simple frequency spectrum.

    Very simple calculation of frequency spectrum with no detrending,
    windowing, etc, just the first half (positive frequency components) of
    abs(fft(x))

    Parameters
    ----------
    x : array_like
        The input array, in the time-domain.

    Returns
    -------
    spec : array_like
        The frequency spectrum of `x`.

    """
    spec = np.absolute(np.fft.fft(x))
    spec = spec[:len(x) / 2]  # take positive frequency components
    spec /= len(x)  # normalize
    spec *= 2.0  # to get amplitudes of sine components, need to multiply by 2
    spec[0] /= 2.0  # except for the dc component
    return spec


class TuningCurve(object):
    """Class to facilitate working with tuning curves."""

    def __init__(self, D=None):
        """
        If `D` is a dict, it is used to give initial values to the tuning curve.
        """
        self._tuning_curves = {}
        self._counts = {}
        if D is not None:
            for k,v in D.items():
                self._tuning_curves[k] = [v]
                self._counts[k] = 1
                self.n = 1
        else:
            self.n = 0

    def add(self, D):
        for k,v  in D.items():
            self._tuning_curves[k].append(v)
            self._counts[k] += 1
        self.n += 1

    def __getitem__(self, i):
        D = {}
        for k,v in self._tuning_curves[k].items():
            D[k] = v[i]
        return D

    def __repr__(self):
        return "TuningCurve: %s" % self._tuning_curves

    def stats(self):
        """Return the mean tuning curve with stderrs."""
        mean = {}
        stderr = {}
        n = self.n
        for k in self._tuning_curves.keys():
            arr = np.array(self._tuning_curves[k])
            mean[k] = arr.mean()
            stderr[k] = arr.std()*n/(n-1)/np.sqrt(n)
        return mean, stderr

    def max(self):
        """Return the key of the max value and the max value."""
        k = _dict_max(self._tuning_curves)
        return k, self._tuning_curves[k]
