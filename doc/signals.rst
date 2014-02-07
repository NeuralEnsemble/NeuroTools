The ``signals`` module
======================

This module provides key functions to manipulate and perform calculations
on data structures like spike trains, spike lists, analog signals.

The import and export of these data structures is typically achieved using the NeuroTools `io` module.

Spike Trains
------------

A spike train is a sorted vector of spike times, which is the result of a simulation or acquired by measurement.
It has therefore some attributes, like ``t_start`` and ``t_stop``, which must in general
be specified by the user, as they can not be inferred from the data.

**Note:** the standard time unit used by NeuroTools is milliseconds.

When several spike trains are gathered, they are collected in a ``SpikeList`` object, which is effectively a dictionary
of spike trains, the id of the cells being used as a key. See the ``SpikeList`` class for more details.


The ``SpikeTrain`` class
------------------------

Creation
~~~~~~~~

Here are some examples of creating simple ``SpikeTrain`` objects::

    >>> spk1 = SpikeTrain(arange(0,100,10), t_start = 0, t_stop=100)
    >>> spk2 = SpikeTrain(arange(0,100,10))

As you can see, defining ``t_start`` and ``t_stop`` is optional. If those attributes are not provided,
they are inferred from the data as the min and the max of the spike times::

    >>> spk1.t_stop
        100
    >>> spk2.t_stop
        90

You can access the raw vector of the spike times with the ``spike_times`` attribute::

    >>> spk1.spike_times
        array([  0.,  10.,  20.,  30.,  40.,  50.,  60.,  70.,  80.,  90.])

The ``SpikeTrain`` object by itself may be interesting, but it will be more often used within a
``SpikeList`` object, as we'll see later.

The ``SpikeTrain`` object has many useful methods such as ``mean_rate()``, ``copy()``, ``isi()``, ``cv_isi()``, ``raster_plot()``. See the ``signals`` API documentation for more
details.

An example of basic use of a ``SpikeTrain`` object is as follows::

    >>> spk1.raster_plot()  # Generates a raster  plot between t_start and t_stop
    >>> spk1.isi()
        array([ 10.,  10.,  10.,  10.,  10.,  10.,  10.,  10.,  10.])
    >>> spk1. mean_rate()
        100.0


.. _SpikeList:

The ``SpikeList`` class
-----------------------

A ``SpikeList`` is effectively a dictionary of ``SpikeTrains``, and can be thought of conceptually as the output of a network simulation,
when one records the spikes of a given population.  The ``SpikeList`` is an object to organize such a recording
as a tuple of key-value pairs ``( spike_source_id, spike_train)``, that is, a ``SpikeList`` object acts as a dictionary created as follows:
``{'id' : SpikeTrain, ...}``

Creation
~~~~~~~~

The constructor of a ``SpikeList`` is works as follows: it will accepts a list of tuples ``(id, spike_time)``,
parameters ``t_start`` and ``t_stop``, and the list of all the recorded ids. The last three parameters
can't be inferred from the data safely, so it's better if they are specified by the user.  Here is an simple example of how to create a ``SpikeList``::

    >>> list = [(numpy.random.random_integers(0,10,1)[0],1000*numpy.random.rand()) for i in xrange(1000)]
    >>> spklist = SpikeList(list, range(11), t_start=0, t_stop=1000)

Here, ``range(11)`` is used to specify that in the lists we will have cells with ids between 0 and 10 (this is needed
because we can have silent cells in the ``SpikeList``), and that we will keep spikes only between ``t_start`` and ``t_stop``.
All the ``SpikeTrain`` objects within the ``SpikeList`` will share the same ``t_start`` and ``t_stop``.

Rather than calling the ``SpikeList`` constructor,  a more common way to create a ``SpikeList`` in NeuroTools is to use the ``load_spikelist()`` or the ``load()`` functions.
If you have generated your data with PyNN_, you can use the loading functions made for this purpose. For example if
you have recorded the spikes of a population in a file "spikes.dat", then one can load it as a ``SpikeList`` as follows::

    >>> spklist = load("spikes.dat",'spikes')

Using this syntax, the header information contained in the file is used to create the population, and ``t_start`` and ``t_stop`` are
inferred automatically as the min and the max of all the ``SpikeTrains`` within the ``SpikeList``.
If you want to keep the control on the parameters while creating the ``SpikeList``, do the following::

    >>> spklist = load_spikelist("spikes.dat", range(11), t_start=0, t_stop=1000)


Now that the ``SpikeList`` has been created, you can start to explore it.

Navigation
~~~~~~~~~~

You can access ``SpikeTrain`` objects within the ``Spikelist`` with the simple syntax ``spklist[id]``::

    >>> spklist[0]
        <NeuroTools.signals.SpikeTrain object at 0x92aac0c>
    >>> spklist[0].mean_rate()
        23.5
    >>> for spiketrain in spklist:
            print spktrain.isi()

As you can see in the example, one can navigate and iterate over a ``SpikeList`` object and have access to
all the ``SpikeTrain``\s within the object. To have an explicit list of all the ids contained in the ``SpikeList``,
use the function ``id_list()``::

    >>> spklist.id_list()
        array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    >>> spklist[15]
        <type 'exceptions.Exception'>: id 15 is not present in the SpikeList. See id_list()

You can't access ``SpikeTrains`` of non recorded cells :-)


Slices
~~~~~~

You can do slices of your ``SpikeList`` object, either according to the time axis or to selected ids::

    >>> subspklist = spklist.time_slice(500,1500)
    >>> subspklist = spklist.id_slice(50) # Will select 50 random id within spklist.id_list()
    >>> subspklist = spklist.id_slice([2,3,5,6])

If you want to select only cells matching a particular criteria within the ``SpikeList``, you can use the
``select_ids`` method::

    >>> subspklist = spklist.select_ids("cell.mean_rate() > 0")
    >>> subspklist = spklist.select_ids("mean(cell.isi()) < 1")



Viewing and saving
~~~~~~~~~~~~~~~~~~

Several signal object methods have the property of being able to generate plots. These plots are implemented to be generic and customizable.  As can be seen in the API documentation, plotting functions accept a ``display`` flag, which you can set to ``True`` if you want to generate a new figure, or you can set it to a subplot if you want to do more complicated plots.  You can also provided additional
parameters to the plot function of pylab as extra arguments::

    >>> spklist.raster_plot(100, t_start=0, t_stop=500)
    >>> spklist.raster_plot(100, t_start=0, t_stop=500, display=subplot(221), kwargs={'color':'red'})

To save the ``SpikeList``, you can use either a Standard Text file output, or a Pickle (compressed) file. To do that, just
provide to the ``save()`` method the corresponding file objects::

    >>> spklist.save("spikes.dat") # Default mode, will create a text file
    >>> spklist.save(StandPickleFile("spikes.dat"))


--------------
Analog Signals
--------------

``NeuroTools.signals`` also handles analog signals. These are generally also
recorded during a simulation or experiment, such as for example a Vm trace, a conductance or a current. Such a signal is
defined by a number of values between ``t_start`` and ``t_stop`` with a time step ``dt``.

The ``AnalogSignal`` class
--------------------------

Creation
~~~~~~~~

When we create an ``AnalogSignal``, we have to provide the list of the data, the time step of their acquisition,
and as an option the ``t_start`` and ``t_stop`` parameters. If ``None``, ``t_start`` will be 0 and ``t_stop`` will be ``len(data)/dt``::

    >>> x = AnalogSignal(sin(arange(1000)),0.1)
    >>> x.t_stop
        100
    >>> x = AnalogSignal(sin(arange(1000)),0.1,0,50)
    >>> x.t_stop
        50
    >>> len(x)
        500

You can access the raw data of the ``AnalogSignal`` by just using::

    >>> x.signal

Usage
~~~~~

Several functions can be applied to an ``AnalogSignal``. See the global API for an exhaustive list. You can for example
do an ``event_trigger_average()``, slice the signal according to some events, detect areas of the signals which are
over a certain threshold, and so on...


The ``AnalogSignalList`` class
------------------------------

As for the ``SpikeList``, the ``AnalogSignalList`` is a collection of ``AnalogSignal`` objects. It has the same structure as
the ``SpikeList``, meaning this is a dictionary containing ``AnalogSignal``\s with the key being the id of the cells.

Creation
~~~~~~~~

The constructor of an ``AnalogSignalList`` is made as follows: it will accept an array with all the signals, the time step,
and additional parameters like ``t_start``, ``t_stop``, and the list of all the recorded ids. The last three parameters
can't be inferred from the data safely, so it's better if they are specified by the user. Nevertheless, the most
common way to create ``AnalogSignal``\s in NeuroTools is to use the ``load_analogsignal()`` or ``load()`` functions, as
explained below. Currently, the constructor of the ``AnalogSignalList`` is mainly tuned to be used with these load functions,
and it is therefore not so simple to create one from a list of ``AnalogSignal``\s::

    >>> sig1 = AnalogSignal(sin(arange(10000),dt=0.1,t_start=0, t_stop=1000)
    ... There seems to be something missing here?


Navigation
~~~~~~~~~~

You can access an ``AnalogSignal`` object within the ``AnalogSignalList`` with the simple syntax ``aslist[id]``::

    >>> aslist[0]
        <NeuroTools.signals.AnalogSignal object at 0x92aac0c>
    >>> mean(aslist[0].signal)
        23.5
    >>> for as in aslist:
            print as.signal

As you can see in the example, one can navigate and iterate over an ``AnalogSignalList`` object and have access to
all the ``AnalogSignal``\s within the object. To have an explicit list of all the ids contained in the ``AnalogSignalList``,
use the function ``id_list()``::

    >>> aslist.id_list()
        array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    >>> aslist[15]
        <type 'exceptions.Exception'>: id 15 is not present in the AnalogSignalList. See id_list()

You can't access ``AnalogSignal``\s of non recorded cells :-)

Viewing and saving
~~~~~~~~~~~~~~~~~~

Similarly to the ``SpikeList object``, some methods of the ``AnalogSignalList`` object can
generate plots. Again, a plot can be either a new figure or a subplot::

    >>> vmlist.plot(1, display= subplot(221), kwargs={'color':'r'})

Again, also, those ``AnalogSignalList`` objects can be saved to file. You just have to use the already implemented ``FileHandler`` (see io.py) ``StandardTextFile`` (default) or ``StandardPickleFile`` and use the ``save()`` method::

    >>> vmlist.save("vm.dat") # Will create a text file by default
    >>> vmlist.save(StandardTextFile("vm.dat")) # save as before
    >>> vmlist.save(StandardPickleFile("vm.dat")) # Saved as pickled object
    >>> vm = load(StandardPickleFile("vm.dat"),'v) # Read the pickle file

Note that for the moment, there is a slight distinction for the conductance files, since the ``load`` function is
tuned for PyNN. Since PyNN saves exc/inh conductances in the same file, the ``load`` function, called on a file
generated by PyNN, will return two ``AnalogSignalList`` ::

    >>> ge, gi = load("conductance.dat",'g')


Autodoc
-------

.. automodule:: NeuroTools.signals
   :members:
   :undoc-members:

    
.. _PyNN: http://neuralensemble.org/PyNN/
