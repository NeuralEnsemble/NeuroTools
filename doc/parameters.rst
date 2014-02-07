The ``parameters`` module
=========================

We consider it to be best practice to cleanly separate the parameters of a model
from the model itself. At the least, parameters should be defined in a separate
section at the start of a file. Ideally, they should be defined in a separate
file entirely. This makes version control easier, since the model code typically
changes less often than the parameters, and makes it easier to track a
simulation project, since the parameter sets can be stored in a database,
displayed in a GUI, etc.


Parameters
----------

At their simplest, individual parameters consist of a name and a value. The
value is either a simple type such as a numerical value or a string, or an
aggregate of such simple types, such as a set, list or array.

However, we may also wish to specify the physical dimensions of the parameter,
i.e., its units, and the range of permissible values.

It is also often useful to specify an object that generates numerical values or
strings, such as a random number generator, and treat that object as the
parameter.

To support all these uses, we define the ``Parameter`` and ``ParameterRange``
classes, and various subclasses of the ``ParameterDist`` abstract class, such as
``GammaDist``, ``NormalDist`` and ``UniformDist``.


The ``Parameter`` class
-----------------------

Here are some examples of creating ``Parameter`` objects::

    >>> i1 = Parameter(3)
    >>> f1 = Parameter(6.2)
    >>> f2 = Parameter(-65.3, "mV")
    >>> s1 = Parameter("hello", name="message_to_the_world")
    
The parameter name, units, value and type can be accessed as attributes::

    >>> i1.value
    3
    >>> f1.type
    <type 'float'>
    >>> f2.units
    'mV'
    >>> s1.name
    'message_to_the_world'

``Parameter`` objects are not hugely useful at the moment. The units are not
used for checking dimensional consistency, for example, and ``Parameter``
objects are not drop-in replacements for numerical values - you must always use
the ``value`` attribute to access the value, whereas it might be nice to define,
for example, a class ``IntegerParameter`` which was a subclass of the built-in
``int`` type.


The ``ParameterRange`` class
----------------------------

When investigating the behaviour of a model or in doing sensitivity analysis, it
is often useful to run a model several times using a different value for a
certain parameter each time (also see the ``iter_range_keys()`` and similar
methods of the ``ParameterSet`` class, below). The ``ParameterRange`` class
supports this. Some usage examples::

    >>> tau_m_range = ParameterRange([10.0, 15.0, 20.0], "ms", "tau_m")
    >>> tau_m_range.name
    'tau_m'
    >>> tau_m_range.next()
    10.0
    >>> tau_m_range.next()
    15.0
    >>> [2*tau_m for tau_m in tau_m_range]
    [20.0, 30.0, 40.0]


The ``ParameterDist`` classes
-----------------------------

As with taking parameter values from a series or range, it is often useful to
pick values from a particular random distribution. Three classes are available:
``UniformDist``, ``GammaDist`` and ``NormalDist``. Examples::
    
    >>> ud = UniformDist(min=-1.0, max=1.0)
    >>> gd = GammaDist(mean=0.5, std=1.0)
    >>> nd = NormalDist(mean=-70, std=5.0)
    >>> ud.next()
    array([-0.56342352]) 
    >>> gd.next(3)
    array([ 0.04061142,  0.05550265,  0.23469344])
    >>> nd.next(2)
    array([-76.18506715, -68.71229944]) 

[Note that very similar functionality is available with the ``RandomDistribution``
class in the `pyNN.random` module. We should look at the best ways to avoid
duplication].

Parameter sets
--------------

A problem with parameter sets for large-scale, detailed models is that the list
of parameters gets very long and unwieldy, and due to the typically hierarchical
nature of such models, the individual parameter names can also get very long,
e.g., ``v1_layer5_pyramidal_apical_dend_gbar_na``.

A solution to this is to give the parameter set a hierarchical structure as well,
which allows the top-level list of parameters to be very short (e.g. ``v1``,
``retina`` and ``lgn`` for a visual system simulation) since the top-level
parameters are themselves parameter sets.

The simplest way to implement this in Python is using nested dicts. One
disadvantage of this is that accessing deeply-nested parameters can be very
verbose, e.g. ``v1['layer5']['pyramidal']['apical_dend']['na']['gbar']``. A
second disadvantage is that it is tedious to flatten the hierarchy when
this becomes necessary, e.g. for serialisation - writing to file, etc.

For these reasons we have created a ``ParameterSet`` class, which:

  1. allows a more convenient notation;
  
  2. enables subsets of the parameters, lower in the hierarchy, to be passed
  around by themselves;
  
  3. provides convenient methods for reading from/writing to file and for
  determining the differences between two different parameter sets.
  
An example of the notation is ``v1.layer5.pyramidal.apical_dend.na.gbar``, which
requires only a single `.` for each level in the hierarchy rather than two
"``'``"s, a "``[``" and a "``]``". This is not much shorter than
``v1_layer5_pyramidal_apical_dend_gbar_na`` - the difference is that
``v1.layer5.pyramidal`` is itself a ``ParameterSet`` object that can be passed
as an argument to the pyramidal cell object, which doesn't care about
``v1.layer4.spinystellate``, let alone ``retina.ganglioncell.magno.tau_m``
(while ``v1_layer5_pyramidal`` is just a ``NameError``).

The ``ParameterSet`` class
--------------------------

Creation
~~~~~~~~

``ParameterSet`` objects may be created from a dict::

    >>> sim_params = ParameterSet({'dt': 0.11, 'tstop': 1000.0})
    
or loaded from a URL::

    >>> exc_cell_params = ParameterSet("https://neuralensemble.org/svn/NeuroTools/trunk/doc/example.param")

They may be nested::

    >>> inh_cell_params = ParameterSet({'tau_m': 15.0, 'cm': 0.5})
    >>> network_params = ParameterSet({'excitatory_cells': exc_cell_params, 'inhibitory_cells': inh_cell_params})
    >>> P = ParameterSet({'sim': sim_params, 'network': network_params}, label="my_params")

Note that although we show here only numerical parameter values,
``Parameter``, ``ParameterRange`` and ``ParameterDist`` objects, as well as
strings, may also be parameter values.

Navigation
~~~~~~~~~~

Individual parameters may be accessed/set using dot notation::

    >>> P.sim.dt
    0.11
    >>> P.network.inhibitory_cells.tau_m
    15.0
    >>> P.network.inhibitory_cells.cm = 0.75
    
or the usual dictionary access notation::

    >>> P['network']['inhibitory_cells']['cm']
    0.75
    
or mixing the two (which may be required if some of the parameter names contain
spaces)::

    >>> P['network'].excitatory_cells['tau_m']
    10.0

Viewing and saving
~~~~~~~~~~~~~~~~~~

To see the entire parameter set at once, nicely formatted use the ``pretty()``
method::

    >>> print P.pretty()
    {
      "network": {
        "excitatory_cells": url("https://neuralensemble.org/svn/NeuroTools/trunk/doc/example.param")
        "inhibitory_cells": {
          "tau_m": 15.0,
          "cm": 0.75,
        },
      },
      "sim": {
        "tstop": 1000.0,
        "dt": 0.11,
      },
    }

By default, if the ``ParameterSet`` contains other ``ParameterSet``\s that were
loaded from URLs, these will be represented with a ``url()`` function in the
output, but there is also the option to expand all URLs and show the full
contents::

    >>> print P.pretty(expand_urls=True)
    {
      "network": {
        "excitatory_cells": {
          "tau_refrac": 0.11,
          "tau_m": 10.0,
          "cm": 0.25,
          "synI": {
            "tau": 10.0,
            "E": -75.0,
          },
          "synE": {
            "tau": 1.5,
            "E": 0.0,
          },
          "v_thresh": -57.0,
          "v_reset": -70.0,
          "v_rest": -70.0,
        },
        "inhibitory_cells": {
          "tau_m": 15.0,
          "cm": 0.75,
        },
      },
      "sim": {
        "tstop": 1000.0,
        "dt": 0.11,
      },
    }

If a ``ParameterSet`` was loaded from a URL, it may be modified then saved back
to the same URL, provided the protocol supports writing::

    >>> exc_cell_params.save()
    Traceback (most recent call last):
      File "<stdin>", line 1, in ?
      File "parameters.py", line 266, in save
        raise Exception("Saving using the %s protocol is not implemented" % scheme)
    Exception: Saving using the https protocol is not implemented
    
or saved to a different URL::

    >>> exc_cell_params.save(url="file:///tmp/exc_params")

The file format is the same as that produced by the ``pretty()`` method.

Copying and converting
~~~~~~~~~~~~~~~~~~~~~~

A ``ParameterSet`` can be used simply as a dictionary, but can also be
converted explicitly to a ``dict`` if required::

    >>> print sim_params.as_dict()
    {'tstop': 1000.0, 'dt': 0.11}

[need to say something about ``tree_copy()``]

Iteration
~~~~~~~~~

There are several different ways to iterate over all or part of the
``ParameterSet`` object. ``keys()``, ``values()`` and ``items()`` work as for
``dict``s. For the sake of more readable code, ``names()`` is provided as an
alias for ``keys()`` and ``parameters()`` as an alias for ``items()``::

    >>> P.names()
    ['network', 'sim']
    >>> exc_cell_params.parameters()
    [('tau_refrac', 0.11), ('tau_m', 10.0), ('cm', 0.25),
     ('synI', {'tau': 10.0, 'E': -75.0}), ('synE', {'tau': 1.5, 'E': 0.0}),
     ('v_thresh', -57.0), ('v_reset', -70.0), ('v_rest', -70.0)]
    
To flatten nested parameter sets, i.e., the iterate recursively over all
branches of the tree, the the ``flatten()`` method returns a ``dict`` with keys
created by joining the names at each hierarchical level with a separator
character ('.' by default)::

    >>> network_params.flatten()
    {'excitatory_cells.synI.E': -75.0, 'excitatory_cells.v_rest': -70.0,
     'excitatory_cells.tau_refrac': 0.11, 'excitatory_cells.v_reset': -70.0,
     'excitatory_cells.v_thresh': -57.0, 'excitatory_cells.tau_m': 10.0,
     'excitatory_cells.synI.tau': 10.0, 'excitatory_cells.cm': 0.25,
     'inhibitory_cells.cm': 0.75, 'excitatory_cells.synE.tau': 1.5,
     'excitatory_cells.synE.E': 0.0, 'inhibitory_cells.tau_m': 15.0}

while the ``flat()`` method returns a generator which yields
``(name, value)`` tuples.::

    >>> for x in network_params.flat():
    ...   print x
    ...


The ``ParameterSpace`` class
----------------------------

The ``ParameterSpace`` class is a subclass of ``ParameterSet`` that is
allowed to contain ``ParameterRange`` and ``ParameterDist`` objects as
parameters. This turns the single point in parameter space represented by a
``ParameterSet`` into a set of points. For example, the following definition
creates a set of six points in parameter space, which can be obtained in turn
using the ``iter_inner()`` method::

    >>> PS = ParameterSpace({
    ...        'x': 999,
    ...        'y': ParameterRange([10, 20]),
    ...        'z': ParameterRange([-1, 0, 1])
    ... })
    >>> for P in PS.iter_inner():
    ...     print P
    {'y': 10, 'x': 999, 'z': -1}
    {'y': 20, 'x': 999, 'z': -1}
    {'y': 10, 'x': 999, 'z': 0}
    {'y': 20, 'x': 999, 'z': 0}
    {'y': 10, 'x': 999, 'z': 1}
    {'y': 20, 'x': 999, 'z': 1}

Putting parameter distribution objects inside a ``ParameterSpace`` allows an
essentially infinite number of points to be generated::

    >>> PS2 = ParameterSpace({
    ...    'x': UniformDist(min=-1.0, max=1.0),
    ...    'y': GammaDist(mean=0.5, std=1.0),
    ...    'z': NormalDist(mean=-70, std=5.0)
    ... })
    >>> for P in PS2.realize_dists(n=3):
    ...     print P
    {'y': 1.81311773668, 'x': 0.883293989399, 'z': -73.5871002759}
    {'y': 0.299391158731, 'x': 0.371474054049, 'z': -68.6936045978}
    {'y': 2.90108202422, 'x': -0.388218831787, 'z': -68.6681724449}


Autodoc
-------

.. automodule:: NeuroTools.parameters
   :members:
   :undoc-members:


