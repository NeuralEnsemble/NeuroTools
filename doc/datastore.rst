The ``datastore`` module
========================

The ``datastore`` module aims to present a consistent interface for persistent
data storage, irrespective of storage back-end. The main intended usage is for
caching of intermediate results. If it takes a long computation time to
calculate an object's data, it can save the data to a datastore the first time
it is run, and then if an identical object is needed in future, it can just
retrieve its data from the store and not have to compute it.

Since it is intended for objects to be able to store part or all of their
internal data, the storage/retrieval keys are based on the object identity and
state.

We assume that an object's identity is uniquely defined by its type (which may
also depend on the source code revision number) and its parameters, while its
state is defined by its identity and by its inputs (we should possibly add some
concept of time to this).

Hence, any object (which we call a 'component' in this context) must have
the following attributes:

``parameters``
  a NeuroTools ``ParameterSet`` object

``input``
  another component or ``None``; we assume a single input.
  A list of inputs should also be possible. We need to be wary of
  recurrent loops, in which two components both have each other as
  direct or indirect inputs).

``version``
  the source-code version

There are two advantages to using the ``datastore`` module rather than just
using, say, ``shelve`` directly::

1. You don't have to worry about keeping track of the key used to identify
   your data in the store: the ``DataStore`` object takes care of this for
   you.
2. You can use various different back-ends to store your data (local
   filesystem, remote filesystem, database) and to manage the keys
   (``shelve``, a database, the filesystem), and the interface remains the
   same.

Creating a datastore
--------------------

Two different storage backends are available, ``ShelveDataStore`` and
``DjangoORMDataStore``. It is also intended to
be easy to write your own, custom storage backend. Whichever backend is used,
after you have created your datastore, the interface is the same. For this
example we will use the ``ShelveDataStore``::

    >>> from NeuroTools.datastore import ShelveDataStore
    >>> datastore = ShelveDataStore(root_dir="/tmp")
 
Here we specify that the ``shelve`` files will be created in ``/tmp``. Now let
us create a simple component whose data we wish to store::

    >>> class SimpleComponent(object):
    ...     def __init__(self, parameters,

Autodoc
-------

.. automodule:: NeuroTools.datastore
   :members:
   :undoc-members:
