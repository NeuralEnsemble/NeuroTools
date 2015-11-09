# -*- coding: utf-8 -*-
"""
NeuroTools.signals.spikes
==================

A collection of functions to create, manipulate and play with spikes signals. 

Classes
-------

SpikeTrain       - object representing a spike train, for one cell. Useful for plots, 
                   calculations such as ISI, CV, mean rate(), ...
SpikeList        - object representing the activity of a population of neurons. Functions as a
                   dictionary of SpikeTrain objects, with methods to compute firing rate,
                   ISI, CV, cross-correlations, and so on.

Functions
---------

load_spikelist       - load a SpikeList object from a file. Expects a particular format.
                       Can also load data in a different format, but then you have
                       to write your own File object that will know how to read the data (see io.py)
load                 - a generic loader for all the previous load methods.

See also NeuroTools.signals.analogs
"""

import os, re, numpy
import scipy.signal
import logging
from NeuroTools import check_dependency, check_numpy_version
from NeuroTools import analysis

from NeuroTools.io import *
from NeuroTools.plotting import get_display, set_axis_limits, set_labels, SimpleMultiplot, progress_bar
from .pairs import *
from .intervals import *


from NeuroTools import check_dependency
HAVE_MATPLOTLIB = check_dependency('matplotlib')
if HAVE_MATPLOTLIB:
    import matplotlib
    matplotlib.use('Agg')
HAVE_PYLAB = check_dependency('pylab')
if HAVE_PYLAB:
    import pylab
else:
    PYLAB_ERROR = "The pylab package was not detected"
if not HAVE_MATPLOTLIB:
    MATPLOTLIB_ERROR = "The matplotlib package was not detected"

newnum = check_numpy_version()

# check whether numpy.histogram supports the deprecated "new" keyword.
hist_new = True
try:
    numpy.histogram(numpy.array([1,2,3]), new=True)
except TypeError:
    hist_new = False

class SpikeTrain(object):
    """
    SpikeTrain(spikes_times, t_start=None, t_stop=None)
    This class defines a spike train as a list of times events.

    Event times are given in a list (sparse representation) in milliseconds.

    Inputs:
        spike_times - a list/numpy array of spike times (in milliseconds)
        t_start     - beginning of the SpikeTrain (if not, this is infered)
        t_stop      - end of the SpikeTrain (if not, this is infered)

    Examples:
        >> s1 = SpikeTrain([0.0, 0.1, 0.2, 0.5])
        >> s1.isi()
            array([ 0.1,  0.1,  0.3])
        >> s1.mean_rate()
            8.0
        >> s1.cv_isi()
            0.565685424949
    """

    #######################################################################
    ## Constructor and key methods to manipulate the SpikeTrain objects  ##
    #######################################################################
    def __init__(self, spike_times, t_start=None, t_stop=None):
        #TODO: add information about sampling rate at time of creation
        
        """
        Constructor of the SpikeTrain object
        
        See also
            SpikeTrain
        """

        self.t_start     = t_start
        self.t_stop      = t_stop
        self.spike_times = numpy.array(spike_times, numpy.float32)

        # If t_start is not None, we resize the spike_train keeping only
        # the spikes with t >= t_start
        if self.t_start is not None:
            self.spike_times = numpy.extract((self.spike_times >= self.t_start), self.spike_times)

        # If t_stop is not None, we resize the spike_train keeping only
        # the spikes with t <= t_stop
        if self.t_stop is not None:
            self.spike_times = numpy.extract((self.spike_times <= self.t_stop), self.spike_times)

        # We sort the spike_times. May be slower, but is necessary by the way for quite a 
        # lot of methods...
        self.spike_times = numpy.sort(self.spike_times, kind="quicksort")
        # Here we deal with the t_start and t_stop values if the SpikeTrain
        # is empty, with only one element or several elements, if we
        # need to guess t_start and t_stop
        # no element : t_start = 0, t_stop = 0.1 
        # 1 element  : t_start = time, t_stop = time + 0.1
        # several    : t_start = min(time), t_stop = max(time)
        
        size = len(self.spike_times)
        if size == 0:
            if self.t_start is None: 
                self.t_start = 0
            if self.t_stop is None:
                self.t_stop  = 0.1
        elif size == 1: # spike list may be empty
            if self.t_start is None:
                self.t_start = self.spike_times[0]
            if self.t_stop is None:
                self.t_stop = self.spike_times[0] + 0.1
        elif size > 1:
            if self.t_start is None:
                self.t_start = numpy.min(self.spike_times)
            if numpy.any(self.spike_times < self.t_start):
                raise ValueError("Spike times must not be less than t_start")
            if self.t_stop is None:
                self.t_stop = numpy.max(self.spike_times)
            if numpy.any(self.spike_times > self.t_stop):
                raise ValueError("Spike times must not be greater than t_stop")

        if self.t_start >= self.t_stop :
            raise Exception("Incompatible time interval : t_start = %s, t_stop = %s" % (self.t_start, self.t_stop))
        if self.t_start < 0:
            raise ValueError("t_start must not be negative")
        if numpy.any(self.spike_times < 0):
            raise ValueError("Spike times must not be negative")

    def __str__(self):
        return str(self.spike_times)

    def __del__(self):
        del self.spike_times

    def __len__(self):
        return len(self.spike_times)
    
    def __getslice__(self, i, j):
        """
        Return a sublist of the spike_times vector of the SpikeTrain
        """
        return self.spike_times[i:j]
    
    def time_parameters(self):
        """
        Return the time parameters of the SpikeTrain (t_start, t_stop)
        """
        return (self.t_start, self.t_stop)
    
    def is_equal(self, spktrain):
        """
        Return True if the SpikeTrain object is equal to one other SpikeTrain, i.e
        if they have same time parameters and same spikes_times
        
        Inputs:
            spktrain - A SpikeTrain object
        
        See also:
            time_parameters()
        """
        test = (self.time_parameters() == spktrain.time_parameters())
        return numpy.all(self.spike_times == spktrain.spike_times) and test
    
    def copy(self):
        """
        Return a copy of the SpikeTrain object
        """
        return SpikeTrain(self.spike_times, self.t_start, self.t_stop)


    def duration(self):
        """
        Return the duration of the SpikeTrain
        """
        return self.t_stop - self.t_start
    
    
    def merge(self, spiketrain, relative=False):
        """
        Add the spike times from a spiketrain to the current SpikeTrain
        
        Inputs:
            spiketrain - The SpikeTrain that should be added
            relative - if True, relative_times() is called on both spiketrains before merging
        
        Examples:
            >> a = SpikeTrain(range(0,100,10),0.1,0,100)
            >> b = SpikeTrain(range(400,500,10),0.1,400,500)
            >> a.merge(b)
            >> a.spike_times
                [   0.,   10.,   20.,   30.,   40.,   50.,   60.,   70.,   80.,
                90.,  400.,  410.,  420.,  430.,  440.,  450.,  460.,  470.,
                480.,  490.]
            >> a.t_stop
                500
        """
        if relative:
            self.relative_times()
            spiketrain.relative_times()
        self.spike_times = numpy.insert(self.spike_times, self.spike_times.searchsorted(spiketrain.spike_times), \
                                        spiketrain.spike_times)
        self.t_start     = min(self.t_start, spiketrain.t_start)
        self.t_stop      = max(self.t_stop, spiketrain.t_stop)

    def format(self, relative=False, quantized=False):
        """
        Return an array with a new representation of the spike times
        
        Inputs:
            relative  - if True, spike times are expressed in a relative
                       time compared to the previsous one
            quantized - a value to divide spike times with before rounding
            
        Examples:
            >> st.spikes_times=[0, 2.1, 3.1, 4.4]
            >> st.format(relative=True)
                [0, 2.1, 1, 1.3]
            >> st.format(quantized=2)
                [0, 1, 2, 2]
        """
        spike_times = self.spike_times.copy()

        if relative and len(spike_times) > 0:
            spike_times[1:] = spike_times[1:] - spike_times[:-1]

        if quantized:
            assert quantized > 0, "quantized must either be False or a positive number"
            #spike_times =  numpy.array([time/self.quantized for time in spike_times],int)
            spike_times = (spike_times/quantized).round().astype('int')

        return spike_times

    def jitter(self,jitter):
        """
        Returns a new SpikeTrain with spiketimes jittered by a normal distribution.

        Inputs:
              jitter - sigma of the normal distribution

        Examples:
              >> st_jittered = st.jitter(2.0)
        """
        
        return SpikeTrain(self.spike_times+jitter*(numpy.random.normal(loc=0.0,scale=1.0,size=self.spike_times.shape[0])),t_start=self.t_start,t_stop=self.t_stop)
        


    #######################################################################
    ## Analysis methods that can be applied to a SpikeTrain object       ##
    #######################################################################
    
    def isi(self):
        """
        Return an array with the inter-spike intervals of the SpikeTrain
        
        Examples:
            >> st.spikes_times=[0, 2.1, 3.1, 4.4]
            >> st.isi()
                [2.1, 1., 1.3]
        
        See also
            cv_isi
        """
        return numpy.diff(self.spike_times)

    def mean_rate(self, t_start=None, t_stop=None):
        """ 
        Returns the mean firing rate between t_start and t_stop, in Hz
        
        Inputs:
            t_start - in ms. If not defined, the one of the SpikeTrain object is used
            t_stop  - in ms. If not defined, the one of the SpikeTrain object is used
        
        Examples:
            >> spk.mean_rate()
                34.2
        """
        if (t_start == None) & (t_stop == None):
            t_start = self.t_start
            t_stop  = self.t_stop
            idx     = self.spike_times
        else:
            if t_start == None: 
                t_start = self.t_start
            else:
                t_start = max(self.t_start, t_start)
            if t_stop == None: 
                t_stop=self.t_stop
            else:
                t_stop = min(self.t_stop, t_stop)
            idx = numpy.where((self.spike_times >= t_start) & (self.spike_times <= t_stop))[0]
        return 1000.*len(idx)/(t_stop-t_start)

    def cv_isi(self):
        """
        Return the coefficient of variation of the isis.
        
        cv_isi is the ratio between the standard deviation and the mean of the ISI
          The irregularity of individual spike trains is measured by the squared
        coefficient of variation of the corresponding inter-spike interval (ISI)
        distribution normalized by the square of its mean.
          In point processes, low values reflect more regular spiking, a
        clock-like pattern yields CV2= 0. On the other hand, CV2 = 1 indicates
        Poisson-type behavior. As a measure for irregularity in the network one
        can use the average irregularity across all neurons.
        
        http://en.wikipedia.org/wiki/Coefficient_of_variation
        
        See also
            isi, cv_kl
            
        """
        isi = self.isi()
        if len(isi) > 0:
            return numpy.std(isi)/numpy.mean(isi)
        else:
            logging.debug("Warning, a CV can't be computed because there are not enough spikes")
            return numpy.nan

    def cv_kl(self, bins=100):
        """
        Provides a measure for the coefficient of variation to describe the
        regularity in spiking networks. It is based on the Kullback-Leibler
        divergence and decribes the difference between a given
        interspike-interval-distribution and an exponential one (representing
        poissonian spike trains) with equal mean.
        It yields 1 for poissonian spike trains and 0 for regular ones.
        
        Reference:
            http://invibe.net/LaurentPerrinet/Publications/Voges08fens
        
        Inputs:
            bins - the number of bins used to gather the ISI
        
        Examples:
            >> spklist.cv_kl(100)
                0.98
        
        See also:
            cv_isi
            
        """
        isi = self.isi() / 1000.
        if len(isi) < 2:
            logging.debug("Warning, a CV can't be computed because there are not enough spikes")
            return numpy.nan
        else:
            if not hist_new:
                proba_isi, xaxis = numpy.histogram(isi, bins=bins, normed=True)
                xaxis = xaxis[:-1]
            else:
                if newnum:
                    proba_isi, xaxis = numpy.histogram(isi, bins=bins, normed=True, new=True)
                    xaxis = xaxis[:-1]
                else:
                    proba_isi, xaxis = numpy.histogram(isi, bins=bins, normed=True)
            proba_isi /= numpy.sum(proba_isi)
            bin_size = xaxis[1]-xaxis[0]
            # differential entropy: http://en.wikipedia.org/wiki/Differential_entropy
            KL = - numpy.sum(proba_isi * numpy.log(proba_isi+1e-16)) + numpy.log(bin_size)
            KL -= -numpy.log(self.mean_rate()) + 1.
            CVkl=numpy.exp(-KL)
            return CVkl
    

    def fano_factor_isi(self):
        """ 
        Return the fano factor of this spike trains ISI.
        
        The Fano Factor is defined as the variance of the isi divided by the mean of the isi
        
        http://en.wikipedia.org/wiki/Fano_factor
        
        See also
            isi, cv_isi
        """
        isi = self.isi()
        if len(isi) > 0:
            fano = numpy.var(isi)/numpy.mean(isi)
            return fano
        else: 
            raise Exception("No spikes in the SpikeTrain !")

    def time_axis(self, time_bin=10):
        """
        Return a time axis between t_start and t_stop according to a time_bin
        
        Inputs:
            time_bin - the bin width
            
        Examples:
            >> st = SpikeTrain(range(100),0.1,0,100)
            >> st.time_axis(10)
                [ 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        
        See also
            time_histogram
        """
        if newnum:
            axis = numpy.arange(self.t_start, self.t_stop+time_bin, time_bin)
        else:
            axis = numpy.arange(self.t_start, self.t_stop, time_bin)
        return axis

    def raster_plot(self, t_start=None, t_stop=None, interval=None, display=True, kwargs={}):
        """
        Generate a raster plot with the SpikeTrain in a subwindow of interest,
        defined by t_start and t_stop. 
        
        Inputs:
            t_start - in ms. If not defined, the one of the SpikeTrain object is used
            t_stop  - in ms. If not defined, the one of the SpikeTrain object is used
            display - if True, a new figure is created. Could also be a subplot
            kwargs  - dictionary contening extra parameters that will be sent to the plot 
                      function
        
        Examples:
            >> z = subplot(221)
            >> st.raster_plot(display=z, kwargs={'color':'r'})
        
        See also
            SpikeList.raster_plot
        """
        if t_start is None: t_start = self.t_start
        if t_stop is None:  t_stop = self.t_stop
        
        if interval is None:
            interval = Interval(t_start, t_stop)
        
        spikes  = interval.slice_times(self.spike_times)
        subplot = get_display(display)
        if not subplot or not HAVE_PYLAB:
            print(PYLAB_ERROR)
        else:
            if len(spikes) > 0:
                subplot.plot(spikes,numpy.ones(len(spikes)),',', **kwargs)
                xlabel = "Time (ms)"
                ylabel = "Neurons #"
                set_labels(subplot, xlabel, ylabel)
                pylab.draw()

    def time_offset(self, offset):
        """
        Add an offset to the SpikeTrain object. t_start and t_stop are
        shifted from offset, so does all the spike times.
         
        Inputs:
            offset - the time offset, in ms
        
        Examples:
            >> spktrain = SpikeTrain(arange(0,100,10))
            >> spktrain.time_offset(50)
            >> spklist.spike_times
                [  50.,   60.,   70.,   80.,   90.,  100.,  110.,  
                120.,  130.,  140.]
        """
        self.t_start += offset
        self.t_stop  += offset
        self.spike_times += offset

    def time_slice(self, t_start, t_stop):
        """ 
        Return a new SpikeTrain obtained by slicing between t_start and t_stop,
        where t_start and t_stop may either be single values or sequences of
        start and stop times.
        
        Inputs:
            t_start - begining of the new SpikeTrain, in ms.
            t_stop  - end of the new SpikeTrain, in ms.

        Examples:
            >> spk = spktrain.time_slice(0,100)
            >> spk.t_start
                0
            >> spk.t_stop
                100
            >>> spk = spktrain.time_slice([20,70], [40,90])
            >>> spk.t_start
                20
            >>> spk.t_stop
                90
            >>> len(spk.time_slice(41, 69))
                0
        """
        if hasattr(t_start, '__len__'):
            if len(t_start) != len(t_stop):
                raise ValueError("t_start has %d values and t_stop %d. They must be of the same length." % (len(t_start), len(t_stop)))
            mask = False
            for t0,t1 in zip(t_start, t_stop):
                mask = mask | ((self.spike_times >= t0) & (self.spike_times <= t1))
            t_start = t_start[0]
            t_stop = t_stop[-1]
        else:
            mask = (self.spike_times >= t_start) & (self.spike_times <= t_stop)
        spikes = numpy.extract(mask, self.spike_times)
        return SpikeTrain(spikes, t_start, t_stop)

    def interval_slice(self, interval):
        """ 
        Return a new SpikeTrain obtained by slicing with an Interval. The new 
        t_start and t_stop values of the returned SpikeTrain are the extrema of the Interval
        
        Inputs:
            interval - The interval from which spikes should be extracted

        Examples:
            >> spk = spktrain.time_slice(0,100)
            >> spk.t_start
                0
            >> spk.t_stop
                100
        """
        times           = interval.slice_times(self.spike_times)
        t_start, t_stop = interval.time_parameters()
        return SpikeTrain(times, t_start, t_stop)
        

    def time_histogram(self, time_bin=10, normalized=True, binary=False):
        """
        Bin the spikes with the specified bin width. The first and last bins
        are calculated from `self.t_start` and `self.t_stop`.
        
        Inputs:
            time_bin   - the bin width for gathering spikes_times
            normalized - if True, the bin values are scaled to represent firing rates
                         in spikes/second, otherwise otherwise it's the number of spikes 
                         per bin.
            binary     - if True, a binary matrix of 0/1 is returned
        
        Examples:
            >> st=SpikeTrain(range(0,100,5),0.1,0,100)
            >> st.time_histogram(10)
                [200, 200, 200, 200, 200, 200, 200, 200, 200, 200]
            >> st.time_histogram(10, normalized=False)
                [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        
        See also
            time_axis
        """
        bins = self.time_axis(time_bin)
        if newnum:
            if not hist_new:
                hist, edges = numpy.histogram(self.spike_times, bins)
            else:
                hist, edges = numpy.histogram(self.spike_times, bins, new=newnum)
        else:
            hist, edges = numpy.histogram(self.spike_times, bins)
        hist = hist.astype(float)
        if normalized: # what about normalization if time_bin is a sequence?
            hist *= 1000.0/float(time_bin)
        if binary:
            hist = hist.astype(bool).astype(int)
        return hist
    
    def instantaneous_rate(self, resolution, kernel, norm, m_idx=None,
                           t_start=None, t_stop=None, acausal=True, trim=False):
        """
        Estimate instantaneous firing rate by kernel convolution.
        
        Inputs:
            resolution  - time stamp resolution of the spike times (ms). the
                          same resolution will be assumed for the kernel
            kernel      - kernel function used to convolve with
            norm        - normalization factor associated with kernel function
                          (see analysis.make_kernel for details)
            t_start     - start time of the interval used to compute the firing
                          rate
            t_stop      - end time of the interval used to compute the firing
                          rate (included)
            acausal     - if True, acausal filtering is used, i.e., the gravity
                          center of the filter function is aligned with the
                          spike to convolve
            m_idx       - index of the value in the kernel function vector that
                          corresponds to its gravity center. this parameter is
                          not mandatory for symmetrical kernels but it is
                          required when assymmetrical kernels are to be aligned
                          at their gravity center with the event times
            trim        - if True, only the 'valid' region of the convolved
                          signal are returned, i.e., the points where there
                          isn't complete overlap between kernel and spike train
                          are discarded
                          NOTE: if True and an assymetrical kernel is provided
                          the output will not be aligned with [t_start, t_stop]
                          
        See also:
            analysis.make_kernel
        """
        
        if t_start is None:
            t_start = self.t_start
            
        if t_stop is None:
            t_stop = self.t_stop
        
        if m_idx is None:
            m_idx = kernel.size / 2
        
        time_vector = numpy.zeros((t_stop - t_start)/resolution + 1)
        
        spikes_slice = self.spike_times[(self.spike_times >= t_start) & (
            self.spike_times <= t_stop)]
        
        for spike in spikes_slice:
            index = (spike - t_start) / resolution
            time_vector[index] = 1
        
        r = norm * scipy.signal.fftconvolve(time_vector, kernel, 'full')

        if acausal is True:
            if trim is False:
                r = r[m_idx:-(kernel.size - m_idx)]
                t_axis = numpy.linspace(t_start, t_stop, r.size)
                return t_axis, r
            
            elif trim is True:
                r = r[2 * m_idx:-2*(kernel.size - m_idx)]
                t_start = t_start + m_idx * resolution
                t_stop = t_stop - ((kernel.size) - m_idx) * resolution
                t_axis = numpy.linspace(t_start, t_stop, r.size)
                return t_axis, r
            
        if acausal is False:
            if trim is False:
                r = r[m_idx:-(kernel.size - m_idx)]
                t_axis = (numpy.linspace(t_start, t_stop, r.size) +
                          m_idx * resolution)
                return t_axis, r
            
            elif trim is True:
                r = r[2 * m_idx:-2*(kernel.size - m_idx)]
                t_start = t_start + m_idx * resolution
                t_stop = t_stop - ((kernel.size) - m_idx) * resolution
                t_axis = (numpy.linspace(t_start, t_stop, r.size) +
                          m_idx * resolution)
                return t_axis, r
    
    def relative_times(self):
        """
        Rescale the spike times to make them relative to t_start.
        
        Note that the SpikeTrain object itself is modified, t_start 
        is subtracted from spike_times, t_start and t_stop
        """
        if self.t_start != 0:
            self.spike_times -= self.t_start
            self.t_stop      -= self.t_start
            self.t_start      = 0.0

    def distance_victorpurpura(self, spktrain, cost=0.5):
        """
        Function to calculate the Victor-Purpura distance between two spike trains.
        See J. D. Victor and K. P. Purpura,
            Nature and precision of temporal coding in visual cortex: a metric-space
            analysis.,
            J Neurophysiol,76(2):1310-1326, 1996
        
        Inputs:
            spktrain - the other SpikeTrain
            cost     - The cost parameter. See the paper for more information
        """
        nspk_1      = len(self)
        nspk_2      = len(spktrain)
        if cost == 0: 
            return abs(nspk_1-nspk_2)
        elif cost > 1e9 :
            return nspk_1+nspk_2
        scr = numpy.zeros((nspk_1+1,nspk_2+1))
        scr[:,0] = numpy.arange(0,nspk_1+1)
        scr[0,:] = numpy.arange(0,nspk_2+1)
            
        if nspk_1 > 0 and nspk_2 > 0:
            for i in xrange(1, nspk_1+1):
                for j in xrange(1, nspk_2+1):
                    scr[i,j] = min(scr[i-1,j]+1,scr[i,j-1]+1)
                    scr[i,j] = min(scr[i,j],scr[i-1,j-1]+cost*abs(self.spike_times[i-1]-spktrain.spike_times[j-1]))
        return scr[nspk_1,nspk_2]


    def distance_kreuz(self, spktrain, dt=0.1):
        """
        Function to calculate the Kreuz/Politi distance between two spike trains
        See  Kreuz, T.; Haas, J.S.; Morelli, A.; Abarbanel, H.D.I. & Politi, A. 
            Measuring spike train synchrony. 
            J Neurosci Methods, 165:151-161, 2007

        Inputs:
            spktrain - the other SpikeTrain
            dt       - the bin width used to discretize the spike times
        
        Examples:
            >> spktrain.KreuzDistance(spktrain2)
        
        See also
            VictorPurpuraDistance
        """
        N              = (self.t_stop-self.t_start)/dt
        vec_1          = numpy.zeros(N, numpy.float32)
        vec_2          = numpy.zeros(N, numpy.float32)
        result         = numpy.zeros(N, float)
        idx_spikes     = numpy.array(self.spike_times/dt,int)
        previous_spike = 0
        if len(idx_spikes) > 0:
            for spike in idx_spikes[1:]:
                vec_1[previous_spike:spike] = (spike-previous_spike)
                previous_spike = spike
        idx_spikes     = numpy.array(spktrain.spike_times/dt,int)
        previous_spike = 0
        if len(idx_spikes) > 0:
            for spike in idx_spikes[1:]:
                vec_2[previous_spike:spike] = (spike-previous_spike)
                previous_spike = spike
        idx = numpy.where(vec_1 < vec_2)[0]
        result[idx] = vec_1[idx]/vec_2[idx] - 1
        idx = numpy.where(vec_1 > vec_2)[0]
        result[idx] = -vec_2[idx]/vec_1[idx] + 1
        return numpy.sum(numpy.abs(result))/len(result)
    
    
    def psth(self, events, time_bin=2, t_min=50, t_max=50, display = False, kwargs={}, average=True):
        """
        Return the psth of the spike times contained in the SpikeTrain according to selected events, 
        on a time window t_spikes - tmin, t_spikes + tmax
        
        Inputs:
            events  - Can be a SpikeTrain object (and events will be the spikes) or just a list 
                      of times
            time_bin- The time bin (in ms) used to gather the spike for the psth
            t_min   - Time (>0) to average the signal before an event, in ms (default 0)
            t_max   - Time (>0) to average the signal after an event, in ms  (default 100)
            display - if True, a new figure is created. Could also be a subplot.
            kwargs  - dictionary contening extra parameters that will be sent to the plot 
                      function
            
        Examples:
            >> spk.psth(spktrain, t_min = 50, t_max = 150)
            >> spk.psth(spktrain, )
            >> spk.psth(range(0,1000,10), display=True)
            
        See also
            SpikeTrain.spike_histogram
        """
        
        if isinstance(events, SpikeTrain):
            events = events.spike_times
        assert (t_min >= 0) and (t_max >= 0), "t_min and t_max should be greater than 0"
        assert len(events) > 0, "events should not be empty and should contained at least one element"

        spk_hist = self.time_histogram(time_bin)
        subplot  = get_display(display)
        count    = 0
        t_min_l  = numpy.floor(t_min/time_bin)
        t_max_l  = numpy.floor(t_max/time_bin)
        result   = numpy.zeros((t_min_l+t_max_l), numpy.float32)
        t_start  = numpy.floor(self.t_start/time_bin)
        t_stop   = numpy.floor(self.t_stop/time_bin)
        result   = []
        for ev in events:
           ev = numpy.floor(ev/time_bin)
           if ((ev - t_min_l )> t_start) and (ev + t_max_l ) < t_stop:
               count  += 1
               result += [spk_hist[(ev-t_min_l):ev+t_max_l]]
        result = numpy.array(result)    
        if average:        
            result /= count
        
        if not subplot or not HAVE_PYLAB:
            return result
        else:
            xlabel = "Time (ms)"
            ylabel = "PSTH"
            time   = numpy.linspace(-t_min, t_max, (t_min+t_max)/time_bin)
            set_labels(subplot, xlabel, ylabel)
            subplot.plot(time, result, c='k', **kwargs)
            xmin, xmax, ymin, ymax = subplot.axis()
            subplot.plot([0,0],[ymin, ymax], c='r')
            set_axis_limits(subplot, -t_min, t_max, ymin, ymax)
            pylab.draw()
        return result
    

    ####################################################################
    ### TOO SPECIFIC METHOD ?
    ### Better documentation
    ####################################################################
    def tuning_curve(self, var_array, normalized=False, method='sum'):
        """
        Calculate a firing-rate tuning curve with respect to some variable.
        Assumes that some variable, such as stimulus orientation, varies
        throughout the recording. The values taken by this variable should be
        supplied in a numpy array `var_array`. The spike train is binned
        according to the number of values in `var_array`, e.g., if there are
        N values, the spikes are binned with a bin width
            (`self.t_stop`-`self.t_start`)/N
        so that each bin is associated with a particular value of the variable
        in `var_array`.

        The return value is a dictionary whose keys are the distinct values in
        `val_array`. The values in the dict depend on the arguments `method` and
        `normalized`.

        If `normalized` is False, the responses (bin values) are spike counts,
        if True, they are firing rates.
        If `method` == "max", each value is the maximum response for a given
        value of the variable.
        If `method` == "sum", each value is the summed response...
        If `method` == "mean", ... you get the idea.

        (If someone can rewrite this more clearly, please do so!)
        """
        binwidth = (self.t_stop - self.t_start)/len(var_array)
        time_histogram = self.time_histogram(binwidth, normalized=normalized)
        assert len(time_histogram) == len(var_array)
        tuning_curve = {}
        counts = {}
        for k, x in zip(var_array, time_histogram):
            if not tuning_curve.has_key(k):
                tuning_curve[k] = 0
                counts[k] = 0
            if method in ('sum', 'mean'):
                tuning_curve[k] += x
                counts[k] += 1
            elif method == 'max':
                tuning_curve[k] = max(x, tuning_curve[k])
            else:
                raise Exception()
        if method == 'mean':
            for k in tuning_curve.keys():
                tuning_curve[k] = tuning_curve[k]/float(counts[k])
        return tuning_curve


    ####################################################################
    ### TOO SPECIFIC METHOD ?
    ### Better documentation
    ####################################################################
    def frequency_spectrum(self, time_bin):
        """
        Returns the frequency spectrum of the time histogram together with the
        frequency axis.
        """
        hist       = self.time_histogram(time_bin, normalized=False)
        from NeuroTools import analysis
        freq_spect = analysis.simple_frequency_spectrum(hist)
        freq_bin   = 1000.0/self.duration() # Hz
        freq_axis  = numpy.arange(len(freq_spect)) * freq_bin    
        return freq_spect, freq_axis    


    ####################################################################
    ### TOO SPECIFIC METHOD ?
    ### Better documentation
    ####################################################################
    def f1f0_ratio(self, time_bin, f_stim):
        """
        Returns the F1/F0 amplitude ratio where the input stimulus frequency is
        f_stim.
        """
        freq_spect = self.frequency_spectrum(time_bin)[0]
        F0 = freq_spect[0]
        freq_bin = 1000.0/self.duration()
        try:
            F1 = freq_spect[int(round(f_stim/freq_bin))]
        except IndexError:
            errmsg = "time_bin (%f) too large to estimate F1/F0 ratio for an input frequency of %f" % (time_bin, f_stim)
            errmsg += "\nFrequency_spectrum: %s" % freq_spect
            raise Exception(errmsg)
        return F1/F0
        


class SpikeList(object):
    """
    SpikeList(spikes, id_list, t_start=None, t_stop=None, dims=None)
    
    Return a SpikeList object which will be a list of SpikeTrain objects.

    Inputs:
        spikes  - a list of (id,time) tuples (id being in id_list)
        id_list - the list of the ids of all recorded cells (needed for silent cells)
        t_start - begining of the SpikeList, in ms. If None, will be infered from the data
        t_stop  - end of the SpikeList, in ms. If None, will be infered from the data
        dims    - dimensions of the recorded population, if not 1D population
    
    t_start and t_stop are shared for all SpikeTrains object within the SpikeList
    
    Examples:
        >> sl = SpikeList([(0, 0.1), (1, 0.1), (0, 0.2)], range(2))
        >> type( sl[0] )
            <type SpikeTrain>
    
    See also
        load_spikelist
    """
    #######################################################################
    ## Constructor and key methods to manipulate the SpikeList objects   ##
    #######################################################################
    def __init__(self, spikes, id_list, t_start=None, t_stop=None, dims=None):
        """
        Constructor of the SpikeList object

        See also
            SpikeList, load_spikelist
        """
        self.t_start     = t_start
        self.t_stop      = t_stop
        self.dimensions  = dims
        self.spiketrains = {}
        id_list          = numpy.sort(id_list)
        
        ##### Implementaion base on pure Numpy arrays, that seems to be faster for
        ## large spike files. Still not very efficient in memory, because we are not
        ## using a generator to build the SpikeList...
        
        if not isinstance(spikes, numpy.ndarray): # is an array:
            spikes = numpy.array(spikes, numpy.float32)
        N = len(spikes)
        
        if N > 0:
            idx    = numpy.argsort(spikes[:,0])
            spikes = spikes[idx]
            logging.debug("sorted spikes[:10,:] = %s" % str(spikes[:10,:]))
        
            break_points = numpy.where(numpy.diff(spikes[:, 0]) > 0)[0] + 1
            break_points = numpy.concatenate(([0], break_points))
            break_points = numpy.concatenate((break_points, [N]))
            for idx in xrange(len(break_points)-1):
                id = spikes[break_points[idx], 0]
                if id in id_list:
                    self.spiketrains[id] = SpikeTrain(spikes[break_points[idx]:break_points[idx+1], 1], self.t_start, self.t_stop)
        
        self.complete(id_list)
        
        if len(self) > 0 and (self.t_start is None or self.t_stop is None):
            self.__calc_startstop()
        
        del spikes

    def __del__(self):
        for id in self.id_list:
            del self.spiketrains[id]

    @property
    def id_list(self):
        """ 
        Return the list of all the cells ids contained in the
        SpikeList object
        
        Examples
            >> spklist.id_list
                [0,1,2,3,....,9999]
        """
        return numpy.array(self.spiketrains.keys(), int)

    def copy(self):
        """
        Return a copy of the SpikeList object
        """
        spklist = SpikeList([], [], self.t_start, self.t_stop, self.dimensions)
        for id in self.id_list:
            spklist.append(id, self.spiketrains[id])
        return spklist

    def __calc_startstop(self):
        """
        t_start and t_stop are shared for all neurons, so we take min and max values respectively.
        TO DO : check the t_start and t_stop parameters for a SpikeList. Is it commun to
        all the spikeTrains within the spikelist or each spikelistes do need its own.
        """
        if len(self) > 0:
            if self.t_start is None:
                start_times = numpy.array([self.spiketrains[idx].t_start for idx in self.id_list], numpy.float32)
                self.t_start = numpy.min(start_times)
                logging.debug("Warning, t_start is infered from the data : %f" %self.t_start)
                for id in self.spiketrains.keys():
                    self.spiketrains[id].t_start = self.t_start
            if self.t_stop is None:
                stop_times = numpy.array([self.spiketrains[idx].t_stop for idx in self.id_list], numpy.float32)
                self.t_stop  = numpy.max(stop_times)
                logging.debug("Warning, t_stop  is infered from the data : %f" %self.t_stop)
                for id in self.spiketrains.keys():
                    self.spiketrains[id].t_stop = self.t_stop
        else:
            raise Exception("No SpikeTrains")

    def __getitem__(self, id):
        if id in self.id_list:
            return self.spiketrains[id]
        else:
            raise Exception("id %d is not present in the SpikeList. See id_list" %id)
    
    def __getslice__(self, i, j):
        """
        Return a new SpikeList object with all the ids between i and j
        """
        ids = numpy.where((self.id_list >= i) & (self.id_list < j))[0]
        return self.id_slice(ids)
    
    #def __setslice__(self, i, j):

    def __setitem__(self, id, spktrain):
        assert isinstance(spktrain, SpikeTrain), "A SpikeList object can only contain SpikeTrain objects"
        self.spiketrains[id] = spktrain
        self.__calc_startstop()

    def __iter__(self):
        return self.spiketrains.itervalues()

    def __len__(self):
        return len(self.spiketrains)

    def __sub_id_list(self, sub_list=None):
        """
        Internal function used to get a sublist for the Spikelist id list
        
        Inputs:
            sublist - can be an int (and then N random cells are selected). Otherwise
                      sub_list is a list of cell in self.id_list. If None, id_list is returned
        
        Examples:
            >> self.__sub_id_list(50)
        """
        if sub_list == None:
            return self.id_list
        elif type(sub_list) == int:
            #import ipdb
            #ipdb.set_trace()
            return numpy.random.permutation(self.id_list)[0:sub_list]
        else:
            return sub_list

            
    def __select_with_pairs__(self, nb_pairs, pairs_generator):
        """
        Internal function used to slice two SpikeList according to a list
        of pairs.  Return a list of pairs
        
        Inputs:
            nb_pairs        - an int specifying the number of cells desired
            pairs_generator - a pairs generator
        
        Examples:
            >> self.__select_with_pairs__(50, RandomPairs(spk1, spk2))
        
        See also
            RandomPairs, AutoPairs, CustomPairs
        """
        pairs  = pairs_generator.get_pairs(nb_pairs)
        spk1   = pairs_generator.spk1.id_slice(pairs[:,0])
        spk2   = pairs_generator.spk2.id_slice(pairs[:,1])
        return spk1, spk2, pairs
   
    def append(self, id, spktrain):
        """
        Add a SpikeTrain object to the SpikeList
        
        Inputs:
            id       - the id of the new cell
            spktrain - the SpikeTrain object representing the new cell
        
        The SpikeTrain object is sliced according to the t_start and t_stop times
        of the SpikeLlist object
        
        Examples
            >> st=SpikeTrain(range(0,100,5),0.1,0,100)
            >> spklist.append(999, st)
                spklist[999]
        
        See also
            concatenate, __setitem__
        """
        assert isinstance(spktrain, SpikeTrain), "A SpikeList object can only contain SpikeTrain objects"
        if id in self.id_list:
            raise Exception("id %d already present in SpikeList. Use __setitem__ (spk[id]=...) instead()" %id)
        else:
            self.spiketrains[id] = spktrain.time_slice(self.t_start, self.t_stop)
            
    def time_parameters(self):
        """
        Return the time parameters of the SpikeList (t_start, t_stop)
        """
        return (self.t_start, self.t_stop)

    def jitter(self,jitter):
        """
        Returns a new SpikeList with spiketimes jittered by a normal distribution.

        Inputs:
              jitter - sigma of the normal distribution

        Examples:
              >> st_jittered = st.jitter(2.0)
        """
        new_SpkList = SpikeList([], [], self.t_start, self.t_stop, self.dimensions)
        for id in self.id_list:
            new_SpkList.append(id, self.spiketrains[id].jitter(jitter))
        return new_SpkList
        

    def time_axis(self, time_bin):
        """
        Return a time axis between t_start and t_stop according to a time_bin
        
        Inputs:
            time_bin - the bin width
        
        See also
            spike_histogram
        """
        if newnum:
            axis = numpy.arange(self.t_start, self.t_stop+time_bin, time_bin)
        else:
            axis = numpy.arange(self.t_start, self.t_stop, time_bin)
        return axis

    def concatenate(self, spklists):
        """
        Concatenation of SpikeLists to the current SpikeList.
        
        Inputs:
            spklists - could be a single SpikeList or a list of SpikeLists
        
        The concatenated SpikeLists must have similar (t_start, t_stop), and
        they can't shared similar cells. All their ids have to be different.
        
        See also
            append, merge, __setitem__
        """
        if isinstance(spklists, SpikeList):
            spklists = [spklists]
        # We check that Spike Lists have similar time_axis
        for sl in spklists:
            if not sl.time_parameters() == self.time_parameters():
                raise Exception("Spike Lists should have similar time_axis")
        for sl in spklists:
            for id in sl.id_list:
                self.append(id, sl.spiketrains[id])

    def merge(self, spikelist, relative=False):
        """
        For each cell id in spikelist that matches an id in this SpikeList,
        merge the two SpikeTrains and save the result in this SpikeList.
        Note that SpikeTrains with ids not in this SpikeList are appended to it.
        
        Inputs:
            spikelist - the SpikeList that should be merged to the current one
            relative  - if True, spike times are expressed in a relative
                        time compared to the previsous one
            
        Examples:
            >> spklist.merge(spklist2)
            
        See also:
            concatenate, append, __setitem__
        """
        for id, spiketrain in spikelist.spiketrains.items():
            if id in self.id_list:
                self.spiketrains[id].merge(spiketrain, relative)
            else:
                if relative:
                    spiketrain.relative_times()
                self.append(id, spiketrain)
        self.__calc_startstop()    
    
    def complete(self, id_list):
        """
        Complete the SpikeList by adding Sempty SpikeTrain for all the ids present in
        ids that will not already be in the SpikeList
        
         Inputs:
            id_list - The id_list that should be completed
        
        Examples:
            >> spklist.id_list
                [0,2,5]
            >> spklist.complete(arange(5))
            >> spklist.id_list
                [0,1,2,3,4]
        """
        id_list     = set(id_list)
        missing_ids = id_list.difference(set(self.id_list))
        for id in missing_ids:
            self.append(id, SpikeTrain([],self.t_start, self.t_stop))
        
    
    def id_slice(self, id_list):
        """
        Return a new SpikeList obtained by selecting particular ids
        
        Inputs:
            id_list - Can be an integer (and then N random cells will be selected)
                      or a sublist of the current ids
        
        The new SpikeList inherits the time parameters (t_start, t_stop)
        
        Examples:
            >> spklist.id_list
                [830, 1959, 1005, 416, 1011, 1240, 729, 59, 1138, 259]
            >> new_spklist = spklist.id_slice(5)
            >> new_spklist.id_list
                [1011, 729, 1138, 416, 59]

        See also
            time_slice, interval_slice
        """
        new_SpkList = SpikeList([], [], self.t_start, self.t_stop, self.dimensions)
        id_list = self.__sub_id_list(id_list)
        for id in id_list:
            try:
                new_SpkList.append(id, self.spiketrains[id])
            except Exception:
                logging.debug("id %d is not in the source SpikeList or already in the new one" %id)
        return new_SpkList

    def time_slice(self, t_start, t_stop):
        """
        Return a new SpikeList obtained by slicing between t_start and t_stop
        
        Inputs:
            t_start - begining of the new SpikeTrain, in ms.
            t_stop  - end of the new SpikeTrain, in ms.
        
        See also
            id_slice, interval_slice
        """
        new_SpkList = SpikeList([], [], t_start, t_stop, self.dimensions)
        for id in self.id_list:
            new_SpkList.append(id, self.spiketrains[id].time_slice(t_start, t_stop))
        new_SpkList.__calc_startstop()
        return new_SpkList
        
    def interval_slice(self, interval):
        """ 
        Return a new SpikeList obtained by slicing with an Interval. The new 
        t_start and t_stop values of the returned SpikeList are the extrema of the Interval
        
        Inputs:
            interval - The Interval to slice with
        
        See also
            id_slice, time_slice
        """
        t_start, t_stop = interval.time_parameters()
        new_SpkList = SpikeList([], [], t_start, t_stop, self.dimensions)
        for id in self.id_list:
            new_SpkList.append(id, self.spiketrains[id].interval_slice(interval))
        return new_SpkList
    
    def time_offset(self, offset):
        """
        Add an offset to the whole SpikeList object. t_start and t_stop are
        shifted from offset, so does all the SpikeTrain.
         
        Inputs:
            offset - the time offset, in ms
        
        Examples:
            >> spklist.t_start
                1000
            >> spklist.time_offset(50)
            >> spklist.t_start
                1050
        """
        self.t_start += offset
        self.t_stop  += offset
        for id in self.id_list:
            self.spiketrains[id].time_offset(offset)
    
    def id_offset(self, offset):
        """
        Add an offset to the whole SpikeList object. All the id are shifted
        according to an offset value.
         
        Inputs:
            offset - the id offset
        
        Examples:
            >> spklist.id_list
                [0,1,2,3,4]
            >> spklist.id_offset(10)
            >> spklist.id_list
                [10,11,12,13,14]
        """
        id_list     = numpy.sort(self.id_list)
        N           = len(id_list)
        if (offset > 0):        
            for idx in xrange(1, len(id_list)+1):
                id  = id_list[N-idx]
                spk = self.spiketrains.pop(id)
                self.spiketrains[id + offset] = spk
        if (offset < 0):        
            for idx in xrange(0, len(id_list)):
                id  = id_list[idx]
                spk = self.spiketrains.pop(id)
                self.spiketrains[id + offset] = spk
    
    def first_spike_time(self):
        """
        Get the time of the first real spike in the SpikeList
        """
        first_spike = self.t_stop
        is_empty    = True
        for id in self.id_list:
            if len(self.spiketrains[id]) > 0:
                is_empty = False
                if self.spiketrains[id].spike_times[0] < first_spike:
                    first_spike = self.spiketrains[id].spike_times[0]
        if is_empty:
            raise Exception("No spikes can be found in the SpikeList object !")
        else:
            return first_spike
    
    def last_spike_time(self):
        """
        Get the time of the last real spike in the SpikeList
        """
        last_spike = self.t_start
        is_empty    = True
        for id in self.id_list:
            if len(self.spiketrains[id]) > 0:
                is_empty = False
                if self.spiketrains[id].spike_times[-1] > last_spike:
                    last_spike = self.spiketrains[id].spike_times[-1]
        if is_empty:
            raise Exception("No spikes can be found in the SpikeList object !")
        else:
            return last_spike
    
    
    def select_ids(self, criteria):
        """
        Return the list of all the cells in the SpikeList that will match the criteria
        expressed with the following syntax. 
        
        Inputs : 
            criteria - a string that can be evaluated on a SpikeTrain object, where the 
                       SpikeTrain should be named ``cell''.
        
        Exemples:
            >> spklist.select_ids("cell.mean_rate() > 0") (all the active cells)
            >> spklist.select_ids("cell.mean_rate() == 0") (all the silent cells)
            >> spklist.select_ids("len(cell.spike_times) > 10")
            >> spklist.select_ids("mean(cell.isi()) < 1")
        """
        selected_ids = []
        for id in self.id_list:
            cell = self.spiketrains[id]
            if eval(criteria):
                selected_ids.append(id)
        return selected_ids

    def sort_by(self, criteria, descending=False):
        """
        Return an array with all the ids of the cells in the SpikeList, 
        sorted according to a particular criteria.
        
        Inputs:
            criteria   - the criteria used to sort the cells. It should be a string
                         that can be evaluated on a SpikeTrain object, where the 
                         SpikeTrain should be named ``cell''.
            descending - if True, then the cells are sorted from max to min. 
            
        Examples:
            >> spk.sort_by('cell.mean_rate()')
            >> spk.sort_by('cell.cv_isi()', descending=True)
            >> spk.sort_by('cell.distance_victorpurpura(target, 0.05)')
        """
        criterias = numpy.zeros(len(self), float)
        for count, id in enumerate(self.id_list):
            cell  = self.spiketrains[id]
            criterias[count] = eval(criteria)
        result = self.id_list[numpy.argsort(criterias)]
        if descending:
            return result[numpy.arange(len(result)-1, -1, -1)]
        else:
            return result

    
    def save(self, user_file):
        """
        Save the SpikeList in a text or binary file
        
        Inputs:
            user_file - The user file that will have its own read/write method
                        By default, if s tring is provided, a StandardTextFile object
                        will be created. Nevertheless, you can also
                        provide a StandardPickleFile
            
        Examples:
            >> spk.save("spikes.txt")
            >> spk.save(StandardTextFile("spikes.txt"))
            >> spk.save(StandardPickleFile("spikes.pck"))
        
        See also:
            DataHandler
        """
        spike_loader = DataHandler(user_file, self)
        spike_loader.save()

    #######################################################################
    ## Analysis methods that can be applied to a SpikeTrain object       ##
    #######################################################################

    
    def isi(self):
        """
        Return the list of all the isi vectors for all the SpikeTrains objects
        within the SpikeList.
        
        See also:
            isi_hist
        """
        isis = []
        for id in self.id_list:
            isis.append(self.spiketrains[id].isi())
        return isis


    def isi_hist(self, bins=50, display=False, kwargs={}):
        """
        Return the histogram of the ISI.
        
        Inputs:
            bins    - the number of bins (between the min and max of the data) 
                      or a list/array containing the lower edges of the bins.
            display - if True, a new figure is created. Could also be a subplot
            kwargs  - dictionary contening extra parameters that will be sent to the plot 
                      function
        
        Examples:
            >> z = subplot(221)
            >> spklist.isi_hist(10, display=z, kwargs={'color':'r'})

        See also:
            isi
        """
        isis          = numpy.concatenate(self.isi())
        if newnum:
            if not hist_new:
                values, xaxis = numpy.histogram(isis, bins=bins)
                xaxis = xaxis[:-1]
            else:
                values, xaxis = numpy.histogram(isis, bins=bins, new=newnum)
                xaxis = xaxis[:-1]
        else:
            values, xaxis = numpy.histogram(isis, bins=bins, new=True)
        subplot       = get_display(display)
        values = values/float(values.sum())
        if not subplot or not HAVE_PYLAB:
            return values, xaxis
        else:
            xlabel = "Inter Spike Interval (ms)"
            ylabel = "Probability"
            set_labels(subplot, xlabel, ylabel)
            subplot.plot(xaxis, values, **kwargs)
            subplot.set_yticks([]) # arbitrary units
            pylab.draw()

        
    def cv_isi(self, float_only=False):
        """
        Return the list of all the CV coefficients for each SpikeTrains object
        within the SpikeList. Return NaN when not enough spikes are present
        
        Inputs:
            float_only - False by default. If true, NaN values are automatically
                         removed
        
        Examples:
            >> spklist.cv_isi()
                [0.2,0.3,Nan,2.5,Nan,1.,2.5]
            >> spklist.cv_isi(True)
                [0.2,0.3,2.5,1.,2.5]

        See also:
            cv_isi_hist, cv_local, cv_kl, SpikeTrain.cv_isi
            
        """
        ids = self.id_list
        N   = len(ids)
        cvs_isi = numpy.empty(N)
        for idx in xrange(N):
            cvs_isi[idx] = self.spiketrains[ids[idx]].cv_isi()

        if float_only:
            cvs_isi = numpy.extract(numpy.logical_not(numpy.isnan(cvs_isi)),cvs_isi)
        return cvs_isi


    def cv_kl(self, bins = 50, float_only=False):
        """
        Return the list of all the CV coefficients for each SpikeTrains object
        within the SpikeList. Return NaN when not enough spikes are present
        
        Inputs:
            bins       - The number of bins used to gathered the ISI
            float_only - False by default. If true, NaN values are automatically
                         removed
        
        Examples:
            >> spklit.cv_kl(50)
                [0.4, Nan, 0.9, nan]
            >> spklist.cv_kl(50, True)
                [0.4, 0.9]

        See also:
            cv_isi_hist, cv_local, cv_isi, SpikeTrain.cv_kl
        """
        ids = self.id_list
        N = len(ids)
        cvs_kl = numpy.empty(N)
        for idx in xrange(N):
            cvs_kl[idx] = self.spiketrains[ids[idx]].cv_kl(bins = bins)
  
        if float_only:
            cvs_kl = numpy.extract(numpy.logical_not(numpy.isnan(cvs_kl)),cvs_kl)
        return cvs_kl

    def cv_isi_hist(self, bins=50, display=False, kwargs={}):
        """
        Return the histogram of the cv coefficients.
        
        Inputs:
            bins    - the number of bins (between the min and max of the data) 
                      or a list/array containing the lower edges of the bins.
            display - if True, a new figure is created. Could also be a subplot
            kwargs  - dictionary contening extra parameters that will be sent to the plot 
                      function
        
        Examples:
            >> z = subplot(221)
            >> spklist.cv_isi_hist(10, display=z, kwargs={'color':'r'})
        
        See also:
            cv_isi, cv_local, cv_kl
        """
        cvs = self.cv_isi(float_only=True)
        if newnum:
            if not hist_new:
                values, xaxis = numpy.histogram(cvs, bins=bins)
                xaxis = xaxis[:-1]
            else:
                values, xaxis = numpy.histogram(cvs, bins=bins, new=newnum)
                xaxis = xaxis[:-1]
        else:
            values, xaxis = numpy.histogram(cvs, bins=bins, new=True)
        subplot       = get_display(display)
        values = values/float(values.sum())                
        if not subplot or not HAVE_PYLAB:
            return values, xaxis
        else:
            xlabel = " CV ISI"
            ylabel = "\% of Neurons"
            set_labels(subplot, xlabel, ylabel)
            subplot.plot(xaxis, values, **kwargs)
            pylab.draw()

    def cv_local(self, t_start=None, t_stop=None, length=12, step=6):
        """
        Provides a modified version of the coefficient of variation, a measure
        that describes the regularity of spiking neurons/networks.
        This CV is local in time, i.e., it considers the isi variance
        for consecutive inter-spike-intervals, and is then averaged
        over a certain window in time (1/12 length), which is shifted in time (1/6 step).
        Therefore rate changes and/or bimodal isi distributions do not lead
        to an exceptionally high CV.
        
        Inputs:
            t_start  - The time to start the averaging 
            t_stop   - The time to stop the averaging 
            length   - A factor to determine the window length for the average, 
                       The window considered will be (t_stop-t_start)/length.
                       12 by default
            time_bin -  factor to determin the step size
        for its shifting.
        
        Examples:
            >> spklist.cv_local(0, 1000, 12, 10)
        
        See also
            cv_isi, cv_isi_hist, cv_kl
        """
        if t_start == None:
            t_start = self.t_start
        if t_stop  == None:
            t_stop  = self.t_stop
        windowLength = (t_stop-t_start)/length
        stepSize = windowLength/step
        maxBin   = int((t_stop-t_start-windowLength)/stepSize)
        #print 'wL=',windowLength,' sSize=',stepSize,' maxBin=',maxBin
        vLocCV   = numpy.zeros(maxBin)
        vCnt     = numpy.zeros(maxBin)
        N = len(self.id_list)
        for i in xrange(N):
            if len(self.spiketrains[i])>15 :
                for j in xrange(2,len(self.spiketrains[i]),1) :
                    diff1=self.spiketrains[i].spike_times[j]-self.spiketrains[i].spike_times[j-1]
                    diff0=self.spiketrains[i].spike_times[j-1]-self.spiketrains[i].spike_times[j-2]
                    tmp=2*numpy.abs(diff1-diff0)/(diff1+diff0)
                    for b in range(maxBin) :
                        if self.spiketrains[i].spike_times[j-2]>b*stepSize+t_start and self.spiketrains[i].spike_times[j]<=b*stepSize+t_start+windowLength :
                            vLocCV[b] = vLocCV[b] + tmp
                            vCnt[b]  = vCnt[b]+1
                        
        locCV= 0.0
        for b in range(maxBin) :
            if vCnt[b] > 0 :
                locCV =  vLocCV[b]/vCnt[b]
    
        return locCV
            

    def mean_rate(self, t_start=None, t_stop=None):
        """
        Return the mean firing rate averaged accross all SpikeTrains between t_start and t_stop.
        
        Inputs:
            t_start - begining of the selected area to compute mean_rate, in ms
            t_stop  - end of the selected area to compute mean_rate, in ms
        
        If t_start or t_stop are not defined, those of the SpikeList are used
        
        Examples:
            >> spklist.mean_rate()
            >> 12.63
        
        See also
            mean_rates, mean_rate_std
        """
        return numpy.mean(self.mean_rates(t_start, t_stop))


    def mean_rate_std(self, t_start=None, t_stop=None):
        """
        Standard deviation of the firing rates accross all SpikeTrains 
        between t_start and t_stop

        Inputs:
            t_start - begining of the selected area to compute std(mean_rate), in ms
            t_stop  - end of the selected area to compute std(mean_rate), in ms
        
        If t_start or t_stop are not defined, those of the SpikeList are used
        
        Examples:
            >> spklist.mean_rate_std()
            >> 13.25

        See also
            mean_rate, mean_rates
        """
        return numpy.std(self.mean_rates(t_start, t_stop))


    def mean_rates(self, t_start=None, t_stop=None):
        """ 
        Returns a vector of the size of id_list giving the mean firing rate for each neuron

        Inputs:
            t_start - begining of the selected area to compute std(mean_rate), in ms
            t_stop  - end of the selected area to compute std(mean_rate), in ms
        
        If t_start or t_stop are not defined, those of the SpikeList are used
        
        See also
            mean_rate, mean_rate_std
        """
        rates = []
        for id in self.id_list:
            rates.append(self.spiketrains[id].mean_rate(t_start, t_stop))
        return rates
    
    
    def rate_distribution(self, nbins=25, normalize=True, display=False, kwargs={}):
        """
        Return a vector with all the mean firing rates for all SpikeTrains.
        
        Inputs:
            bins    - the number of bins (between the min and max of the data) 
                      or a list/array containing the lower edges of the bins.
            display - if True, a new figure is created. Could also be a subplot
            kwargs  - dictionary contening extra parameters that will be sent to the plot 
                      function
        
        See also
            mean_rate, mean_rates
        """
        rates   = self.mean_rates()
        subplot = get_display(display)
        if not subplot or not HAVE_PYLAB:
            return rates
        else:
            if newnum:
                if not hist_new:
                    values, xaxis = numpy.histogram(rates, nbins)
                    xaxis = xaxis[:-1]
                else:
                    values, xaxis = numpy.histogram(rates, nbins, new=newnum)
                    xaxis = xaxis[:-1]
            else:
                values, xaxis = numpy.histogram(rates, nbins)
            xlabel = "Average Firing Rate (Hz)"
            ylabel = "\% of Neurons"
            set_labels(subplot, xlabel, ylabel)
            subplot.plot(xaxis, values/float(values.sum()), **kwargs)
            pylab.draw()
            


    def spike_histogram(self, time_bin, normalized=False, binary=False, display=False, kwargs={}):
        """
        Generate an array with all the spike_histograms of all the SpikeTrains
        objects within the SpikeList.
        
        Inputs:
            time_bin   - the time bin used to gather the data
            normalized - if True, the histogram are in Hz (spikes/second), otherwise they are
                         in spikes/bin
            display    - if True, a new figure is created. Could also be a subplot. The averaged
                         spike_histogram over the whole population is then plotted
            binary     - if True, a binary matrix of 0/1 is returned
            kwargs     - dictionary contening extra parameters that will be sent to the plot 
                         function
        
        See also
            firing_rate, time_axis
        """
        nbins      = self.time_axis(time_bin)
        N          = len(self)
        M          = len(nbins)
        if newnum:
            M -= 1
        if binary:
            spike_hist = numpy.zeros((N, M), numpy.int)
        else:
            spike_hist = numpy.zeros((N, M), numpy.float32)
        subplot    = get_display(display)

        for idx,id in enumerate(self.id_list):
            if newnum:
                if not hist_new:
                    hist, edges = numpy.histogram(self.spiketrains[id].spike_times, nbins)
                else:
                    hist, edges = numpy.histogram(self.spiketrains[id].spike_times, nbins, new=newnum)
            else:
                hist, edges = numpy.histogram(self.spike_times, bins)
            hist = hist.astype(float)
            if normalized: # what about normalization if time_bin is a sequence?
                hist *= 1000.0/float(time_bin)
            if binary:
                hist = hist.astype(bool)
            spike_hist[idx,:] = hist
        if not subplot or not HAVE_PYLAB:
            return spike_hist
        else:
            if normalized:
                ylabel = "Firing rate (Hz)"
            else:
                ylabel = "Spikes per bin"
            xlabel = "Time (ms)"
            set_labels(subplot, xlabel, ylabel)
            axis = self.time_axis(time_bin)
            if newnum:
                axis = axis[:len(axis)-1]
            subplot.plot(axis,numpy.mean(spike_hist, axis=0),**kwargs)
            pylab.draw()
            

    def firing_rate(self, time_bin, display=False, average=False, binary=False, kwargs={}):
        """
        Generate an array with all the instantaneous firing rates along time (in Hz) 
        of all the SpikeTrains objects within the SpikeList. If average is True, it gives the
        average firing rate over the whole SpikeList
        
        Inputs:
            time_bin   - the time bin used to gather the data
            average    - If True, return a single vector of the average firing rate over the whole SpikeList
            display    - if True, a new figure is created. Could also be a subplot. The averaged
                         spike_histogram over the whole population is then plotted
            binary     - If True, a binary matrix with 0/1 is returned. 
            kwargs     - dictionary contening extra parameters that will be sent to the plot 
                         function
        
        See also
            spike_histogram, time_axis
        """
        result = self.spike_histogram(time_bin, normalized=True, binary=False, display=display, kwargs=kwargs)
        if average:
            return numpy.mean(result, axis=0)
        else:
            return result

    def averaged_instantaneous_rate(self, resolution, kernel, norm, m_idx=None,
                                    t_start=None, t_stop=None, acausal=True,
                                    trim=False):
        """
        Estimate the instantaneous firing rate averaged across neurons in the
        SpikeList, by kernel convolution.
        
        Inputs:
            resolution  - time stamp resolution of the spike times (ms). the
                          same resolution will be assumed for the kernel
            kernel      - kernel function used to convolve with
            norm        - normalization factor associated with kernel function
                          (see analysis.make_kernel for details)
            t_start     - start time of the interval used to compute the firing
                          rate
            t_stop      - end time of the interval used to compute the firing
                          rate (included)
            acausal     - if True, acausal filtering is used, i.e., the gravity
                          center of the filter function is aligned with the
                          spike to convolve
            m_idx       - index of the value in the kernel function vector that
                          corresponds to its gravity center. this parameter is
                          not mandatory for symmetrical kernels but it is
                          required when assymmetrical kernels are to be aligned
                          at their gravity center with the event times
            trim        - if True, only the 'valid' region of the convolved
                          signal are returned, i.e., the points where there
                          isn't complete overlap between kernel and spike train
                          are discarded
                          NOTE: if True and an assymetrical kernel is provided
                          the output will not be aligned with [t_start, t_stop]
                          
        See also:
            analysis.make_kernel, SpikeTrain.instantaneous_rate
        """
        
        if t_start is None:
            t_start = self.t_start
            
        if t_stop is None:
            t_stop = self.t_stop
        
        if m_idx is None:
            m_idx = kernel.size / 2
            
        spikes_slice = []
        for i in self:
            train_slice = i.spike_times[(i.spike_times >= t_start) & (
            i.spike_times <= t_stop)]
            spikes_slice += train_slice.tolist()
        
        time_vector = numpy.zeros((t_stop - t_start)/resolution + 1)
        
        for spike in spikes_slice:
            index = (spike - t_start) / resolution
            time_vector[index] += 1.
            
        avg_time_vector = time_vector / float(self.id_list.size)
        
        r = norm * scipy.signal.fftconvolve(avg_time_vector, kernel, 'full')

        if acausal is True:
            if trim is False:
                r = r[m_idx:-(kernel.size - m_idx)]
                t_axis = numpy.linspace(t_start, t_stop, r.size)
                return t_axis, r
            
            elif trim is True:
                r = r[2 * m_idx:-2*(kernel.size - m_idx)]
                t_start = t_start + m_idx * resolution
                t_stop = t_stop - ((kernel.size) - m_idx) * resolution
                t_axis = numpy.linspace(t_start, t_stop, r.size)
                return t_axis, r
            
        if acausal is False:
            if trim is False:
                r = r[m_idx:-(kernel.size - m_idx)]
                t_axis = (numpy.linspace(t_start, t_stop, r.size) +
                          m_idx * resolution)
                return t_axis, r
            
            elif trim is True:
                r = r[2 * m_idx:-2*(kernel.size - m_idx)]
                t_start = t_start + m_idx * resolution
                t_stop = t_stop - ((kernel.size) - m_idx) * resolution
                t_axis = (numpy.linspace(t_start, t_stop, r.size) +
                          m_idx * resolution)
                return t_axis, r
    
    def fano_factor(self, time_bin):
        """
        Compute the Fano Factor of the population activity.
        
        Inputs:
            time_bin   - the number of bins (between the min and max of the data) 
                         or a list/array containing the lower edges of the bins.
        
        The Fano Factor is computed as the variance of the averaged activity divided by its
        mean
        
        See also
            spike_histogram, firing_rate
        """
        firing_rate = self.spike_histogram(time_bin)
        firing_rate = numpy.mean(firing_rate, axis=0)
        fano        = numpy.var(firing_rate)/numpy.mean(firing_rate)
        return fano


    def fano_factors_isi(self):
        """ 
        Return a list containing the fano factors for each neuron
        
        See also
            isi, isi_cv
        """
        fano_factors = []
        for id in self.id_list:
            try:
                fano_factors.append(self.spiketrains[id].fano_factor_isi())
            except:
                pass

        return fano_factors

    
    def id2position(self, id, offset=0):
        """
        Return a position (x,y) from an id if the cells are aranged on a
        grid of size dims, as defined in the dims attribute of the SpikeList object.
        This assumes that cells are ordered from left to right, top to bottom,
        and that dims specifies (height, width), i.e. if dims = (10,12), this is
        an array with 10 rows and 12 columns, and hence has width 12 units and
        height 10 units.
        
        Inputs:
            id - the id of the cell
        
        The 'dimensions' attribute of the SpikeList must be defined
        
        See also
            activity_map, activity_movie
        """
        if self.dimensions is None:
            raise Exception("Dimensions of the population are not defined ! Set spikelist.dimensions")
        if len(self.dimensions) == 1:
            return id-offset
        if len(self.dimensions) == 2:
            x = (id-offset) % self.dimensions[0]
            y = ((id-offset)/self.dimensions[0]).astype(int)
            return (x,y)


    def position2id(self, position, offset=0):
        """
        Return the id of the cell at position (x,y) if the cells are aranged on a
        grid of size dims, as defined in the dims attribute of the SpikeList object
        
        Inputs:
            position - a tuple with the position of the cell
                    
        The 'dimensions' attribute of the SpikeList must be defined and have the same shape
        as the position argument
        
        See also
            activity_map, activity_movie, id2position
        """
        if self.dimensions is None:
            raise Exception("Dimensions of the population are not defined ! Set spikelist.dimensions")
        assert len(position) == len(tuple(self.dimensions)), "position does not have the correct shape !"
        if len(self.dimensions) == 1:
            return position + offset
        if len(self.dimensions) == 2:
            return position[1]*self.dimensions[1] + position[0] + offset


    
    def activity_map(self, t_start=None, t_stop=None, float_positions=None, display=False, kwargs={}):
        """
        Generate a 2D map of the activity averaged between t_start and t_stop.
        If t_start and t_stop are not defined, we used those of the SpikeList object
        
        Inputs:
            t_start         - if not defined, the one of the SpikeList is used
            t_stop          - if not defined, the one of the SpikeList is used
            float_positions - None by default, meaning that the dimensions attribute 
                              of the SpikeList is used to arange the ids on a 2D grid. 
                              Otherwise, if the cells have floating positions, 
                              float_positions should be an array of size
                              (2, nb_cells) with the x (first line) and y (second line) 
                              coordinates of the cells
            display         - if True, a new figure is created. Could also be a subplot. 
                              The averaged spike_histogram over the whole population is 
                              then plotted
            kwargs          - dictionary contening extra parameters that will be sent 
                              to the plot function

        The 'dimensions' attribute of the SpikeList is used to turn ids into 2d positions. It should
        therefore be not empty.
        
        Examples:
            >> spklist.activity_map(0,1000,display=True)
        
        See also
            activity_movie
        """
        subplot = get_display(display)
        
        if t_start == None: 
            t_start = self.t_start
        if t_stop  == None: 
            t_stop  = self.t_stop
        if t_start != self.t_start or t_stop != self.t_stop:
            spklist = self.time_slice(t_start, t_stop)
        else:
            spklist = self
        
        if float_positions is None:
            if self.dimensions is None:
                raise Exception("Dimensions of the population are not defined ! Set spikelist.dims")
            activity_map = numpy.zeros(self.dimensions, float)
            rates        = spklist.mean_rates()
            #id_offset    = min(self.id_list)
            #x,y          = spklist.id2position(spklist.id_list, id_offset)
            x,y          = spklist.id2position(spklist.id_list)            
            #j,i = x, self.dimensions[0] - 1 - y
            for count, id in enumerate(spklist.id_list):
                #activity_map[i[count],j[count]] = rates[count]                            
                activity_map[x[count],y[count]] = rates[count]
            if not subplot or not HAVE_PYLAB or not HAVE_MATPLOTLIB:
                return activity_map
            else:
                im = subplot.imshow(activity_map, **kwargs)
                pylab.colorbar(im)
                pylab.draw()
        elif isinstance(float_positions, numpy.ndarray):
            if not len(spklist.id_list) == len(float_positions[0]):
                raise Exception("Error, the number of flotting positions does not match the number of cells in the SpikeList")
            rates = spklist.mean_rates()
            if not subplot or not HAVE_PYLAB or not HAVE_MATPLOTLIB:
                return rates
            else:
                x = float_positions[0,:]
                y = float_positions[1,:]
                im = subplot.scatter(x,y,c=rates, **kwargs)
                pylab.colorbar(im)
                pylab.draw()


    def pairwise_cc(self, nb_pairs, pairs_generator=None, time_bin=1., average=True, display=False, kwargs={}):
        """
        Function to generate an array of cross correlations computed
        between pairs of cells within the SpikeTrains.
        
        Inputs:
            nb_pairs        - int specifying the number of pairs
            pairs_generator - The generator that will be used to draw the pairs. If None, a default one is
                              created as RandomPairs(spk, spk, no_silent=False, no_auto=True)
            time_bin        - The time bin used to gather the spikes
            average         - If true, only the averaged CC among all the pairs is returned (less memory needed)
            display         - if True, a new figure is created. Could also be a subplot. The averaged
                              spike_histogram over the whole population is then plotted
            kwargs          - dictionary contening extra parameters that will be sent to the plot 
                              function
        
        Examples
            >> a.pairwise_cc(500, time_bin=1, averaged=True)
            >> a.pairwise_cc(500, time_bin=1, averaged=True, display=subplot(221), kwargs={'color':'r'})
            >> a.pairwise_cc(100, CustomPairs(a,a,[(i,i+1) for i in xrange(100)]), time_bin=5)
        
        See also
            pairwise_pearson_corrcoeff, pairwise_cc_zero, RandomPairs, AutoPairs, CustomPairs
        """
        subplot = get_display(display)
        
        ## We have to extract only the non silent cells, to avoid problems
        if pairs_generator is None:
            pairs_generator = RandomPairs(self, self, False, True)

        # Then we select the pairs of cells
        pairs  = pairs_generator.get_pairs(nb_pairs)
        N      = len(pairs)
        if newnum:
            length = 2*(len(pairs_generator.spk1.time_axis(time_bin))-1)
        else:
            length = 2*len(pairs_generator.spk1.time_axis(time_bin))
        if not average:
            results = numpy.zeros((N,length), float)
        else:
            results = numpy.zeros(length, float)
        for idx in xrange(N):
            # We need to avoid empty spike histogram, otherwise the ccf function
            # will give a nan vector
            hist_1 = pairs_generator.spk1[pairs[idx,0]].time_histogram(time_bin)
            hist_2 = pairs_generator.spk2[pairs[idx,1]].time_histogram(time_bin)
            if not average:
                from NeuroTools import analysis
                results[idx,:] = analysis.ccf(hist_1,hist_2)
            else:
                from NeuroTools import analysis
                results += analysis.ccf(hist_1,hist_2)
        if not subplot or not HAVE_PYLAB:
            if not average:
                return results
            else:
                return results/N
        else:
            if average:
                results = results/N
            else:
                results = numpy.sum(results, axis=0)/N
            xaxis   = time_bin*numpy.arange(-len(results)/2, len(results)/2)
            xlabel  = "Time (ms)"
            ylabel  = "Cross Correlation"
            subplot.plot(xaxis, results, **kwargs)
            set_labels(subplot, xlabel, ylabel)
            pylab.draw()


    def pairwise_cc_zero(self, nb_pairs, pairs_generator=None, time_bin=1., time_window=None, display=False, kwargs={}):
        """
        Function to return the normalized cross correlation coefficient at zero time
        lag according to the method given in:
        See A. Aertsen et al, 
            Dynamics of neuronal firing correlation: modulation of effective connectivity
            J Neurophysiol, 61:900-917, 1989
        
        The coefficient is averaged over N pairs of cells. If time window is specified, compute
        the corr coeff on a sliding window, and therefore returns not a value but a vector.
        
        Inputs:
            nb_pairs        - int specifying the number of pairs
            pairs_generator - The generator that will be used to draw the pairs. If None, a default one is
                              created as RandomPairs(spk, spk, no_silent=False, no_auto=True)
            time_bin        - The time bin used to gather the spikes
            time_window     - None by default, and then a single number, the normalized CC is returned.
                              If this is a float, then size (in ms) of the sliding window used to 
                              compute the normalized cc. A Vector is then returned
            display         - if True, a new figure is created. Could also be a subplot. The averaged
                              spike_histogram over the whole population is then plotted
            kwargs          - dictionary contening extra parameters that will be sent to the plot 
                              function
        
        Examples:
            >> a.pairwise_cc_zero(100, time_bin=1)
                1.0
            >> a.pairwise_cc_zero(100, CustomPairs(a, a, [(i,i+1) for i in xrange(100)]), time_bin=1)
                0.45
            >> a.pairwise_cc_zero(100, RandomPairs(a, a, no_silent=True), time_bin=5, time_window=10, display=True)
        
        See also:
            pairwise_cc, pairwise_pearson_corrcoeff, RandomPairs, AutoPairs, CustomPairs
        """
        
        subplot = get_display(display)
        
        if pairs_generator is None:
            pairs_generator = RandomPairs(self, self, False, True)
        
        spk1, spk2, pairs = self.__select_with_pairs__(nb_pairs, pairs_generator)
        N = len(pairs)
        
        if spk1.time_parameters() != spk2.time_parameters():
            raise Exception("The two SpikeList must have common time axis !")
        
        num_bins     = int(numpy.round((self.t_stop-self.t_start)/time_bin)+1)
        mat_neur1    = numpy.zeros((num_bins,N),int)
        mat_neur2    = numpy.zeros((num_bins,N),int)
        times1, ids1 = spk1.convert("times, ids")
        times2, ids2 = spk2.convert("times, ids")
        
        cells_id     = spk1.id_list
        for idx in xrange(len(cells_id)):
            ids1[numpy.where(ids1 == cells_id[idx])[0]] = idx
        cells_id     = spk2.id_list
        for idx in xrange(len(cells_id)):
            ids2[numpy.where(ids2 == cells_id[idx])[0]] = idx
        times1  = numpy.array(((times1 - self.t_start)/time_bin),int)
        times2  = numpy.array(((times2 - self.t_start)/time_bin),int)
        mat_neur1[times1,ids1] = 1
        mat_neur2[times2,ids2] = 1
        if time_window:
            nb_pts   = int(time_window/time_bin)
            mat_prod = mat_neur1*mat_neur2
            cc_time  = numpy.zeros((num_bins-nb_pts),float)
            xaxis    = numpy.zeros((num_bins-nb_pts))
            M        = float(nb_pts*N)
            bound    = int(numpy.ceil(nb_pts/2))
            for idx in xrange(bound,num_bins-bound):
                s_min = idx-bound
                s_max = idx+bound
                Z = numpy.sum(numpy.sum(mat_prod[s_min:s_max]))/M
                X = numpy.sum(numpy.sum(mat_neur1[s_min:s_max]))/M
                Y = numpy.sum(numpy.sum(mat_neur2[s_min:s_max]))/M
                cc_time[s_min] = (Z-X*Y)/numpy.sqrt((X*(1-X))*(Y*(1-Y)))
                xaxis[s_min]   = time_bin*idx
            if not subplot or not HAVE_PYLAB:
                return cc_time
            else:
                xlabel  = "Time (ms)"
                ylabel  = "Normalized CC"
                subplot.plot(xaxis+self.t_start, cc_time, **kwargs)
                set_labels(subplot, xlabel, ylabel)
                pylab.draw()
        else:
            M = float(num_bins*N)
            X = len(times1)/M
            Y = len(times2)/M
            Z = numpy.sum(numpy.sum(mat_neur1*mat_neur2))/M
            return (Z-X*Y)/numpy.sqrt((X*(1-X))*(Y*(1-Y)))


    def distance_victorpurpura(self, nb_pairs, pairs_generator=None, cost=0.5):
        """
        Function to calculate the Victor-Purpura distance averaged over N pairs in the SpikeList
        See J. D. Victor and K. P. Purpura,
            Nature and precision of temporal coding in visual cortex: a metric-space
            analysis.,
            J Neurophysiol,76(2):1310-1326, 1996
        
        Inputs:
            nb_pairs        - int specifying the number of pairs
            pairs_generator - The generator that will be used to draw the pairs. If None, a default one is
                              created as RandomPairs(spk, spk, no_silent=False, no_auto=True)
            cost            - The cost parameter. See the paper for more informations. BY default, set to 0.5
        
        See also
            RandomPairs, AutoPairs, CustomPairs
        """
        if pairs_generator is None:
            pairs_generator = RandomPairs(self, self, False, True)

        pairs = pairs_generator.get_pairs(nb_pairs)
        N     = len(pairs)
        distance   = 0.
        for idx in xrange(N):
            idx_1 = pairs[idx,0]
            idx_2 = pairs[idx,1]
            distance += pairs_generator.spk1[idx_1].distance_victorpurpura(pairs_generator.spk2[idx_2], cost)
        return distance/N
    
    
    def distance_kreuz(self, nb_pairs, pairs_generator=None, dt=0.1):
        """
        Function to calculate the Kreuz/Politi distance between two spike trains
        See Kreuz, T.; Haas, J.S.; Morelli, A.; Abarbanel, H.D.I. & Politi, 
        A. Measuring spike train synchrony. 
        J Neurosci Methods, 2007, 165, 151-161

        Inputs:
            nb_pairs        - int specifying the number of pairs
            pairs_generator - The generator that will be used to draw the pairs. If None, a default one is
                              created as RandomPairs(spk, spk, no_silent=False, no_auto=True)
            dt              - The time bin used to discretized the spike times
        
        See also
            RandomPairs, AutoPairs, CustomPairs
        """
        if pairs_generator is None:
            pairs_generator = RandomPairs(self, self, False, True)

        pairs = pairs_generator.get_pairs(nb_pairs)
        N     = len(pairs)
        
        distance   = 0.
        for idx in xrange(N):
            idx_1 = pairs[idx,0]
            idx_2 = pairs[idx,1]
            distance += pairs_generator.spk1[idx_1].distance_kreuz(pairs_generator.spk2[idx_2], dt)
        return distance/N

    def mean_rate_variance(self, time_bin):
        """
        Return the standard deviation of the firing rate along time,
        if events are binned with a time bin.
        
        Inputs:
            time_bin - time bin to bin events
        
        See also
            mean_rate, mean_rates, mean_rate_covariance, firing_rate
        """
        firing_rate = self.firing_rate(time_bin)
        return numpy.var(numpy.mean(firing_rate, axis=0))

    def mean_rate_covariance(self, spikelist, time_bin):
        """
        Return the covariance of the firing rate along time,
        if events are binned with a time bin.
        
        Inputs:
            spikelist - the other spikelist to compute covariance 
            time_bin  - time bin to bin events
        
        See also
            mean_rate, mean_rates, mean_rate_variance, firing_rate
        """
        if not isinstance(spikelist, SpikeList):
            raise Exception("Error, argument should be a SpikeList object")
        if not spikelist.time_parameters() == self.time_parameters():
            raise Exception("Error, both SpikeLists should share common t_start, t_stop")
        frate_1 = self.firing_rate(time_bin, average=True)
        frate_2 = spikelist.firing_rate(time_bin, average=True)
        N       = len(frate_1)
        cov = numpy.sum(frate_1*frate_2)/N-numpy.sum(frate_1)*numpy.sum(frate_2)/(N*N)
        return cov


    def raster_plot(self, id_list=None, t_start=None, t_stop=None, display=True, kwargs={}):
        """
        Generate a raster plot for the SpikeList in a subwindow of interest,
        defined by id_list, t_start and t_stop. 
        
        Inputs:
            id_list - can be a integer (and then N cells are randomly selected) or a list of ids. If None, 
                      we use all the ids of the SpikeList
            t_start - in ms. If not defined, the one of the SpikeList object is used
            t_stop  - in ms. If not defined, the one of the SpikeList object is used
            display - if True, a new figure is created. Could also be a subplot
            kwargs  - dictionary contening extra parameters that will be sent to the plot 
                      function
        
        Examples:
            >> z = subplot(221)
            >> spikelist.raster_plot(display=z, kwargs={'color':'r'})
        
        See also
            SpikeTrain.raster_plot
        """
        subplot = get_display(display)
        if id_list == None: 
            id_list = self.id_list
            spk = self
        else:
            spk = self.id_slice(id_list)

        if t_start is None: t_start = spk.t_start
        if t_stop is None:  t_stop  = spk.t_stop
        if t_start != spk.t_start or t_stop != spk.t_stop:
            spk = spk.time_slice(t_start, t_stop)

        if not subplot or not HAVE_PYLAB:
            print(PYLAB_ERROR)
        else:
            ids, spike_times = spk.convert(format="[ids, times]")
            idx = numpy.where((spike_times >= t_start) & (spike_times <= t_stop))[0]
            if len(spike_times) > 0:
                subplot.plot(spike_times, ids, ',', **kwargs)
            xlabel = "Time (ms)"
            ylabel = "Neuron #"
            set_labels(subplot, xlabel, ylabel)
            min_id = numpy.min(spk.id_list)
            max_id = numpy.max(spk.id_list)
            length = t_stop - t_start
            set_axis_limits(subplot, t_start-0.05*length, t_stop+0.05*length, min_id-2, max_id+2)
            pylab.draw()


    def psth(self, events, average=True, time_bin=2, t_min=50, t_max=50, display = False, kwargs={}):
        """
        Return the psth of the cells contained in the SpikeList according to selected events, 
        on a time window t_spikes - tmin, t_spikes + tmax
        Can return either the averaged psth (average = True), or an array of all the
        psth triggered by all the spikes.
            
        Inputs:
            events  - Can be a SpikeTrain object (and events will be the spikes) or just a list 
                      of times
            average - If True, return a single vector of the averaged waveform. If False, 
                      return an array of all the waveforms.
            time_bin- The time bin (in ms) used to gather the spike for the psth
            t_min   - Time (>0) to average the signal before an event, in ms (default 0)
            t_max   - Time (>0) to average the signal after an event, in ms  (default 100)
            display - if True, a new figure is created. Could also be a subplot.
            kwargs  - dictionary contening extra parameters that will be sent to the plot 
                      function
            
        Examples:
            >> vm.psth(spktrain, average=False, t_min = 50, t_max = 150)
            >> vm.psth(spktrain, average=True)
            >> vm.psth(range(0,1000,10), average=False, display=True)
            
        See also
            SpikeTrain.spike_histogram
        """
        
        if isinstance(events, SpikeTrain):
            events = events.spike_times
        assert (t_min >= 0) and (t_max >= 0), "t_min and t_max should be greater than 0"
        assert len(events) > 0, "events should not be empty and should contained at least one element"

        spk_hist = self.spike_histogram(time_bin)
        subplot  = get_display(display)
        count    = 0
        t_min_l  = numpy.floor(t_min/time_bin)
        t_max_l  = numpy.floor(t_max/time_bin)
        result   = numpy.zeros((len(self), t_min_l+t_max_l), numpy.float32)
        t_start  = numpy.floor(self.t_start/time_bin)
        t_stop   = numpy.floor(self.t_stop/time_bin)

        for ev in events:
           ev = numpy.floor(ev/time_bin)
           if ((ev - t_min_l )> t_start) and (ev + t_max_l ) < t_stop:               
               count  += 1
               result += spk_hist[:,(ev-t_start-t_min_l):ev-t_start+t_max_l]
        result /= count
            
        if not subplot or not HAVE_PYLAB:
            return result
        else:
            xlabel = "Time (ms)"
            ylabel = "PSTH"
            time   = numpy.linspace(-t_min, t_max, (t_min+t_max)/time_bin)
            set_labels(subplot, xlabel, ylabel)
            if average:
                subplot.plot(time, result, **kwargs)
                subplot.errorbar(times, mean(result, 0), yerr=std(result, 0))
            else:
                for idx in xrange(len(result)):
                    subplot.plot(time, result[idx,:], c='0.5', **kwargs)
                    subplot.hold(1)
                result = numpy.mean(result, 0)
                subplot.plot(time, result, c='k', **kwargs)
            xmin, xmax, ymin, ymax = subplot.axis()
            subplot.plot([0,0],[ymin, ymax], c='r')
            set_axis_limits(subplot, -t_min, t_max, ymin, ymax)
            pylab.draw()
        if average:
            result = numpy.mean(result, 0)
        return result



    def activity_movie(self, time_bin=10, t_start=None, t_stop=None, float_positions=None, output="animation.mpg", bounds=(0,5), fps=10, display=True, kwargs={}):
        """
        Generate a movie of the activity between t_start and t_stop.
        If t_start and t_stop are not defined, we used those of the SpikeList object
        
        Inputs:
            time_bin        - time step to bin activity during the movie. 
                              One frame is the mean rate during time_bin
            t_start         - if not defined, the one of the SpikeList is used, in ms
            t_stop          - if not defined, the one of the SpikeList is used, in ms
            float_positions - None by default, meaning that the dimensions attribute of the SpikeList
                              is used to arange the ids on a 2D grid. Otherwise, if the cells have
                              flotting positions, float_positions should be an array of size
                              (2, nb_cells) with the x (first line) and y (second line) coordinates of
                              the cells
            output          - The filename to store the movie
            bounds          - The common color bounds used during all the movies frame. 
                              This is a tuple
                              of values (min, max), in spikes per frame.
            fps             - The number of frame per second in the final movie
            display         - if True, a new figure is created. Could also be a subplot.
            kwargs          - dictionary contening extra parameters that will be sent to the plot 
                              function

        The 'dimensions' attribute of the SpikeList is used to turn ids into 2d positions. It should
        therefore be not empty.
        
        Examples:
            >> spklist.activity_movie(10,0,1000,bounds=(0,5),display=subplot(221),output="test.mpg")
        
        See also
            activity_map
        """
        subplot = get_display(display)
        if t_start is None: t_start = self.t_start
        if t_stop is None:  t_stop  = self.t_stop
        if not subplot or not HAVE_PYLAB:
            print(PYLAB_ERROR)
        else:
            files        = []
            if float_positions is None:
                activity_map = numpy.zeros(self.dimensions)
                im           = subplot.imshow(activity_map, **kwargs)
                im.set_clim(bounds[0],bounds[1])
                pylab.colorbar(im)
            else:
                rates        = [0]*len(self)
                im           = subplot.scatter(float_positions[0,:], float_positions[1,:], c=rates, **kwargs)
                im.set_clim(bounds[0],bounds[1])
                pylab.colorbar(im)
            count     = 0
            idx       = 0
            manager   = pylab.get_current_fig_manager()
            if t_start != self.t_start or t_stop != self.t_stop:
                spk   = self.time_slice(t_start, t_stop)
            else:
                spk   = self
            time, pos = spk.convert("times, ids")
            # We sort the spikes to allow faster process later
            sort_idx  = time.ravel().argsort(kind="quicksort")
            time      = time[sort_idx]
            pos       = pos[sort_idx]
            x,y       = spk.id2position(pos)
            max_idx   = len(time)-1
            logging.info('Making movie %s - this make take a while' % output)
            if float_positions is None:
                if self.dimensions is None:
                    raise Exception("Dimensions of the population are not defined ! Set spikelist.dims")
                while (t_start < t_stop):
                    activity_map = numpy.zeros(spk.dimensions)
                    while ((time[idx] < t_start + time_bin) and (idx < max_idx)):                                                
                        #j,i = x, self.dimensions[0] - 1 -y
                        activity_map[x[idx],y[idx]] += 1
                        idx += 1
                    im.set_array(activity_map)
                    subplot.title("time = %d ms" %t_start)
                    im.set_clim(bounds[0],bounds[1])
                    manager.canvas.draw()                    
                    fname = "_tmp_spikes_%05d.png" %count
                    #logging.debug("Saving Frame %s", fname)
                    #progress_bar(float(t_start)/t_stop)
                    pylab.savefig(fname)
                    files.append(fname)
                    t_start += time_bin
                    count += 1   
            elif isinstance(float_positions, numpy.ndarray):
                if not len(self) == len(float_positions[0]):
                    raise Exception("Error, the number of flotting positions does not match the number of cells in the SpikeList")
                while (t_start < t_stop):
                    rates = [0]*len(self)
                    while ((time[idx] < t_start + time_bin) and (idx < max_idx)):                                                
                        rates[pos[idx]] += 1
                        idx += 1
                    im = subplot.scatter(float_positions[0,:], float_positions[1,:], c=rates, **kwargs)
                    subplot.title("time = %d ms" %t_start)
                    im.set_clim(bounds[0],bounds[1])
                    manager.canvas.draw()
                    fname = "_tmp_spikes_%05d.png" %count
                    #logging.debug("Saving Frame %s", fname)
                    progress_bar(float(t_start)/t_stop)
                    pylab.savefig(fname)
                    files.append(fname)
                    t_start += time_bin
                    count += 1   
            command = "mencoder 'mf://_tmp_*.png' -mf type=png:fps=%d -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s" %(fps,output)
            logging.debug(command)
            os.system(command)
            ## cleanup
            logging.debug("Clean up....")
            for fname in files: os.remove(fname)


    def pairwise_pearson_corrcoeff(self, nb_pairs, pairs_generator=None, time_bin=1., all_coef=False): 
        """
        Function to return the mean and the variance of the pearson correlation coefficient. 
        For more details, see Kumar et al, ....
        
        Inputs:
            nb_pairs        - int specifying the number of pairs
            pairs_generator - The generator that will be used to draw the pairs. If None, a default one is
                              created as RandomPairs(spk, spk, no_silent=False, no_auto=True)
            time_bin        - The time bin used to gather the spikes
            all_coef        - If True, the whole list of correlation coefficient is returned

        Examples
            >> spk.pairwise_pearson_corrcoeff(50, time_bin=5)
                (0.234, 0.0087)
            >> spk.pairwise_pearson_corrcoeff(100, AutoPairs(spk, spk))
                (1.0, 0.0)
        
        See also
            pairwise_cc, pairwise_cc_zero, RandomPairs, AutoPairs, CustomPairs
        """

        ## We have to extract only the non silent cells, to avoid problems
        if pairs_generator == None:
            pairs_generator = RandomPairs(self, self, False, True)

        pairs = pairs_generator.get_pairs(nb_pairs)
        N     = len(pairs)
        cor   = numpy.zeros(N, float)
        
        for idx in xrange(N):
            hist_1 = pairs_generator.spk1[pairs[idx,0]].time_histogram(time_bin)
            hist_2 = pairs_generator.spk2[pairs[idx,1]].time_histogram(time_bin)
            
            # TODO: normalize the cor, look in 1.6 in Kumar et al.
            # bruederle: the function corrcoeff actually implements the definition in Kumar 1.6
            cov = numpy.corrcoef(hist_1,hist_2)[1][0]

            # bruederle: the expression 'cov' is already, per definition, the pearson correlation coefficient 
            # see http://en.wikipedia.org/wiki/Correlation#The_sample_correlation
            cor[idx] = cov  
            # these two versions have been existing here before
            #cor[count] = cov/numpy.sqrt(n1_hist[0].var()*n2_hist[0].var())

        if all_coef:
            return cor
        else:
            cor_coef_mean = cor.mean()
            cor_coef_std  = cor.std()
        return (cor_coef_mean, cor_coef_std)


    ####################################################################
    ### TOO SPECIFIC METHOD ?
    ### Better documentation
    ####################################################################
    def f1f0_ratios(self, time_bin, f_stim):
        """
        Returns the F1/F0 amplitude ratios for the spike trains contained in the
        spike list, where the input stimulus frequency is f_stim.
        """
        f1f0_dict = {}
        for id, spiketrain in self.spiketrains.items():
            f1f0_dict[id] = spiketrain.f1f0_ratio(time_bin, f_stim)
        return f1f0_dict



    #######################################################################
    ## Method to convert the SpikeList into several others format        ##
    #######################################################################

    def convert(self, format="[times, ids]", relative=False, quantized=False):
        """
        Return a new representation of the SpikeList object, in a user designed format.
            format is an expression containing either the keywords times and ids, 
            time and id.
       
        Inputs:
            relative -  a boolean to say if a relative representation of the spikes 
                        times compared to t_start is needed
            quantized - a boolean to round the spikes_times.
        
        Examples: 
            >> spk.convert("[times, ids]") will return a list of two elements, the 
                first one being the array of all the spikes, the second the array of all the
                corresponding ids
            >> spk.convert("[(time,id)]") will return a list of tuples (time, id)
        
        See also
            SpikeTrain.format
        """
        is_times = re.compile("times")
        is_ids   = re.compile("ids")
        if len(self) > 0:
            times  = numpy.concatenate([st.format(relative, quantized) for st in self.spiketrains.itervalues()])
            ids    = numpy.concatenate([id*numpy.ones(len(st.spike_times), int) for id,st in self.spiketrains.iteritems()])
        else:
            times = []
            ids   = []
        if is_times.search(format):
            if is_ids.search(format):
                return eval(format)
            else:
                raise Exception("You must have a format with [times, ids] or [time, id]")
        is_times = re.compile("time")
        is_ids   = re.compile("id")
        if is_times.search(format):
            if is_ids.search(format):
                result = []
                for id, time in zip(ids, times):
                    result.append(eval(format))
            else:
                raise Exception("You must have a format with [times, ids] or [time, id]")
            return result


    def raw_data(self):
        """
        Function to return a N by 2 array of all times and ids.
        
        Examples:
            >> spklist.raw_data()
            >> array([[  1.00000000e+00,   1.00000000e+00],
                      [  1.00000000e+00,   1.00000000e+00],
                      [  2.00000000e+00,   2.00000000e+00],
                         ...,
                      [  2.71530000e+03,   2.76210000e+03]])
        
        See also:
            convert()
        """
        data = numpy.array(self.convert("[times, ids]"), numpy.float32)
        data = numpy.transpose(data)
        return data



#############################################################
## Object Loaders. Functions used to create NeuroTools
## objects from data generated by pyNN (the most simple form
## supported right now)
#############################################################

def load_spikelist(user_file, id_list=None, t_start=None, t_stop=None, dims=None):
    """
    Returns a SpikeList object from a file. If the file has been generated by PyNN, 
    a header should be found with following parameters:
     ---> dims, dt, id of the first cell, id of the last cell. 
    They must be specified otherwise.  Then the classical PyNN format for text file is:
     ---> one line per event:  absolute time in ms, GID
    
    Inputs:
        user_file - the user_file object with read/write methods. By defaults, if a string
                    is provided, a StandardTextFile object is created
        id_list  - the list of the recorded ids. Can be an int (meaning cells in 
                   the range (0,..,N)), or a list. 
        dims     - if the cells were aranged on a 2/3D grid, a tuple with the dimensions
        t_start  - begining of the simulation, in ms.
        t_stop   - end of the simulation, in ms

    If dims, t_start, t_stop or id_list are None, they will be infered from either 
    the data or from the header. All times are in milliseconds. 
    The format of the file (text, pickle) will be inferred automatically
    """
    spike_loader = DataHandler(user_file)
    return spike_loader.load_spikes(id_list=id_list, t_start=t_start, t_stop=t_stop, dims=dims)


def load(user_file, datatype):
    """
    Convenient data loader for results produced by pyNN. Return the corresponding
    NeuroTools object. Datatype argument may become optionnal in the future, but
    for now it is necessary to specify the type of the recorded data. To have a better control
    on the parameters of the NeuroTools objects, see the load_*** functions.
    
    Inputs:
        user_file - the user_file object with read/write methods. By defaults, if a string
                    is provided, a StandardTextFile object is created
        datatype - A string to specify the type od the data in
                    's' : spikes
                    'g' : conductances
                    'v' : membrane traces
                    'c' : currents
    
    Examples:
        >> load("simulation.dat",'v')
        >> load("spikes.dat",'s')
        >> load(StandardPickleFile("simulation.dat"), 'g')
        >> load(StandardTextFile("test.dat"), 's')
    
    See also:
        load_spikelist, load_conductancelist, load_vmlist, load_currentlist
    """
    if datatype == 's':
        return load_spikelist(user_file)
    elif datatype == 'v':
        return load_vmlist(user_file)
    elif datatype == 'c':
        return load_currentlist(user_file)
    elif datatype == 'g':
        return load_conductancelist(user_file)
    else:
        raise Exception("The datatype %s is not handled ! Should be 's','g','c' or 'v'" %datatype)

def _test():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    #from spikes import *
    _test()
