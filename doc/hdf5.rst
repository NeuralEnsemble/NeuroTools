HDF Tools
=========


HDF5 Overview
-------------

- A `quick start  <http://hdfgroup.com/whatishdf5.html>`_  / `more complete  <http://hdfgroup.com/HDF5/RD100-2002/All_About_HDF5.pdf>`_
-  Hierarchical Data Format version 5, known as HDF5, is a general purpose file format that can store every kind of data in a hierarchical manner. An HDF5 file could be described by a tree of nodes. Each node could contain other nodes or leaves and described with attributes. Leaves under a specified node corresponds to stored data and could also have several attributes to describe them. This format is well adapted to scientific data because it allows to store images, tables, matrices and more with efficient I/O and compression options. That's why it has been chosen as a common file format to share data between NeuroTools developpers. `more details <http://www.pytables.org/docs/manual/>`_

HDF5 and Python
---------------

In fact HDF5 is not only a file format. It's also a high level C++ API that allows to store data in a such format. But it is not the only way to store data in a HDF5 file. The pytables python library, makes easier the conception of HDF5 files in an Object Oriented (OO) manner. It enables files nodes tree navigation, nodes attributes affectation, tables or arrays creation and data insertion. This library is really mature and gives access to many practical functions, but it was necessary to extend this library with new classes and functions to adapt it to Optical Imaging and Neural Network Simulation contexts.

Why developping an extension ?
------------------------------

Developping an extension is an OO approach necessary to make itself reusable, easily extendable and avoid the pytables code pollution.  Moreover, it allows to generate new objects classes that are easier to recognize among base objects and gives access to new functionnalities adapted to each object. A such extension requires that users and developpers who want to work on it have a knowledge on pytables base library. But a class (cf next paragraph) has been fortunately developped to simplify access to pytables classes and function. It will replace old hdf5tools and hdf5api in the NeuroTools folder.

Extension Overview
------------------

This extension is divided in 5 Python files :

#. FileExtension? : the base class of the extension, inherits from the pytables File class. Could be considered as a high level library that gives access to all pytables objects and the new extension objects.
#. AdvancedTable? : add new functionnalities to the pytables Table class like sorts, insertions and queries on tables.
#. Images : useful to convert common file format to VLArray pytables class and display images embedded into HDF5.
#. Movie : encapsulate several common file format images located in a filesystem folder into a VLArray and display them.
#. Spikes, divided into 3 distincts classes :

   #. SpikeList? : store spikelist like [(reltime1, id1), ... (reltimeM, idM)]
   #. SpikingNeurons? : store spikelist like [[reltime1, id1, ..., idM],...,[reltimeK, ..., id N]]
   #. NeuronIndex? : store neuron index defined like [[id1, posX1, posY1, posZ1], ..., [idN, posXN, posYN, posZN]] or [[id1, posX1, posY1, layer1], ..., [idN, posXN, posYN, layerM]]

