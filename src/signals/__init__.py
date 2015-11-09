"""
NeuroTools.signals
==================

A collection of functions to create, manipulate and play with analog signals. 

Classes
-------

AnalogSignal     - object representing an analog signal, with its data. Can be used to do 
                   threshold detection, event triggered averages, ...
AnalogSignalList - list of AnalogSignal objects, again with methods such as mean, std, plot, 
                   and so on
VmList           - AnalogSignalList object used for Vm traces
ConductanceList  - AnalogSignalList object used for conductance traces
CurrentList      - AnalogSignalList object used for current traces

Functions
---------

load_vmlist          - function to load a VmList object (inherits from AnalogSignalList) from a file.
                       Same comments on format as previously.
load_currentlist     - function to load a CurrentList object (inherits from AnalogSignalList) from a file.
                       Same comments on format as previously.
load_conductancelist - function to load a ConductanceList object (inherits from AnalogSignalList) from a file.
                       Same comments on format as previously. load_conductancelist returns two 
                       ConductanceLists, one for the excitatory conductance and one for the inhibitory conductance
load                 - a generic loader for all the previous load methods.
"""

from .spikes import *
from .analogs import *

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
    if datatype in ('s', 'spikes'):
        return load_spikelist(user_file)
    elif datatype == 'v':
        return load_vmlist(user_file)
    elif datatype == 'c':
        return load_currentlist(user_file)
    elif datatype == 'g':
        return load_conductancelist(user_file)
    else:
        raise Exception("The datatype %s is not handled ! Should be 's','g','c' or 'v'" %datatype)
    