The ``io`` module
=================

This module will be the gateway of all the input/output relations in
NeuroTools, especially regarding the inferface with pyNN. This is in
that module that you'll have the Standard Formats currently supported
by NeuroTools (text and pickle, hdf5 planned in a near future), and if
you want to implement your own ``load`` function, reading your own
particular data structure for the ``signals`` module, you should read
the documentation

File Handlers
-------------

A File handler is an abstract object that will have to implement some
key methods in order to be able to read and write NeuroTools objects
from a file (given in the constructor). The idea is that is you want
to design your own File handler, you just have to implement the
abstract methods of the objects, i.e ``write()`` (to write an object
to a file), ``read_spikes(params)`` read data and return a SpikeList
object and ``read_analogs(params, type)``, read data and returns an
analog signal according to type. To have a better understanding, just
have a look to the two file handlers implemented in NeuroTools, i.e
``StandardTextFile`` and ``StandPickleFile``.

The ``StandardTextFile`` class
------------------------------

Creation
~~~~~~~~

The ``StandardTextFile`` inherits from ``FileHandler``

Here is an example of creating simple ``StandardTextFile`` objects::

    >>> textfile = StandardTextFile("test.txt")

Usage
~~~~~

If you want to read a data file with spikes, and return a SpikeList object::
    
    >>> spklist = textfile.read_spikes({'id_list' :range(11), 't_start' : 0, 't_stop' : 1000})

More generally, the ``read_spikes()`` method of an object inheriting from ``FileHandler`` accepts arguments
like id_list, t_start, t_stop, which are the one used in the SpikeList constructor. Note that the ``StandardTextFile`` object have private functions for an internal use only that will check/read 
informations in the headers of the text file, ... See io.py for a deeper understanding of its behavior.

Similar syntax is used for reading a analog signal object::
    
    >>> aslist = textfile.read_analogs('vm', {'id_list':range(11)})

In the case of an ``AnalogSignal``, the type here, selected in [vm, conductance, current] will specified
the type of the NeuroTools object returned by the function. Either a ``VmList``, ``ConductanceList`` or
``CurrentList``

It you want to save an object to a file, just do::
    
    >>> textfile.write(object)

objet can be a SpikeList or any kind of AnalogSignalList.


The ``StandardPickleFile`` class
--------------------------------

Creation
~~~~~~~~

The ``StandardPickleFile`` also inherits from ``FileHandler``

Here is an example of creating simple ``StandardPickleFile`` objects::

    >>> pickfile = StandardPickleFile("test.pick")

Usage
~~~~~

If you want to read a data file with spikes, and return a SpikeList object::
    
    >>> spklist = pickfile.read_spikes({'id_list' : range(11), 't_start' : 0, 't_stop' : 1000})

Since this object inherits from ``FileHandler``, the idea is that its behavior is *exactly* the
same than the ``StandardTextFile``. Similar syntax is used for reading a analog signal object::
    
    >>> aslist = pickfile.read_analogs('vm', {'id_list' : range(11)})

In the case of an ``AnalogSignal``, the type here, selected in [vm, conductance, current] will specified
the type of the NeuroTools object returned by the function. Either a ``VmList``, ``ConductanceList`` or
``CurrentList``

It you want to save an object to a file, just do::
    
    >>> pickfile.write(object)

objet can be a SpikeList or any kind of AnalogSignalList.

The ``YOURStandardFormatFile`` class
------------------------------------

As said before, you just have to implement some key functions, as defined in the ``FileHandler``::

    >>> class YOURStandardFormatFile(FileHandler):
            def write(self, object):
                ### Your method here #########
                ### Should save an object to the file self.filename###
    
            def read_spikes(self, params):
                ### Your method here, reading data from self.filename #########
                ### Should read data and return a SpikeList object constrained by params
                from NeuroTools import signals
                return signals.SpikeList(...)
    
            def read_analogs(self, type, params):
                if not type in ["vm", "current", "conductance"]:
                    raise Exception("The type %s is not available for the Analogs Signals" %type)
                ### Your method here reading data from self.filename #########
                from NeuroTools import signals
                if type == 'vm':
                    return signals.VmList(...)
                elif type == 'conductance':
                    return signals.ConductanceList(...)
                elif type == 'current':
                    return signals.CurrentList(...)

Data Handlers
-------------

The data handler is just a file input/output manager. This is just an interface for ``load/save`` functions.
This is this kind of object which is created by all the ``load`` methods of NeuroTools.signals

The ``DataHandler`` class
-------------------------

You should not have to deal directly with this class, because this is just an interface. See io.py for more details

Autodoc
-------

.. automodule:: NeuroTools.io
   :members:
   :undoc-members:

