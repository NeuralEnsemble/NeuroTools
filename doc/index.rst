NeuroTools
==========

NeuroTools is a collection of tools to support all tasks associated with a neural simulation project and which are not handled by the simulation engine. NeuroTools is written in Python, and works best with `PyNN <http://neuralensemble.org/PyNN>`_, or one of the growing list of simulation engines with a Python front-end such as `NEURON <http://www.neuron.yale.edu/neuron>`_, `NEST <http://www.nest-initiative.org/>`_, `PCSIM <http://sourceforge.net/projects/pcsim/>`_, `FACETS Neuromorphic VLSI <http://wwwasic.kip.uni-heidelberg.de/asicnew/vision/projects/facets/software_and_experiments/>`_, `Brian <http://brian.di.ens.fr/>`_, `MOOSE/GENESIS <http://moose.sourceforge.net/>`_, `Neurospaces/GENESIS <http://neurospaces.sourceforge.net/>`_. NeuroTools provides modules to facilitate simulation setup, parameterization, data management, analysis and visualization. The data-related tools are equally suited to analysis of experimental data, although that is not the primary motivation for their development.

NeuroTools aims to:

#. increase the productivity of individual modellers by automating, simplifying, and establishing best-practices for common tasks,
#. increase the productivity of the neuroscience modelling community by reducing the amount of code duplication across simulation communities,
#. increase the reliability of data analysis tools leveraging `Linus's law <http://en.wikipedia.org/wiki/Given_enough_eyeballs>`_: "given enough eyeballs, all bugs are shallow."

NeuroTools is open-source software, and anyone who is interested is welcome to contribute (`more <http://neuralensemble.org/trac/NeuroTools#Developerpages>`_).

Users' Guide
------------

NeuroTools functionality is modularized as follows:

-  The ``signals`` module provides core functionality that allows manipulation  of and calculations with spike trains and analog signals.
-  The ``spike2`` module offers an easy way for reading data from CED's Spike2  Son files into the NeuroTools enviroment.
-  The ``parameters`` module contains classes to make managing large,  hierarchical parameter sets easier.
-  The ``analysis`` module contains miscellaneous analysis functions
-  The ``stgen`` module contains various stochastic process generators relevant   for Neuroscience (OU, poisson, inhomogenous gamma, ...).
-  The ``utilities`` sub-package containsi the ``srblib``, an easy-to-use  interface for interacting with `SRB   <http://www.sdsc.edu/srb/index.php/Main_Page>`_ servers.
-  The ``io`` module is the gateway for all reading/writing of files, in different formats, in NeuroTools.
-  The ``plotting`` module contains a collection of tools for plotting and image processing, based on Matplotlib and the Python Imaging Library.
-  The ``datastore`` module presents a consistent interface for persistent data storage (e.g., for caching intermediate results), irrespective of storage back-end.

See also more real-world examples in the ``examples`` folder and the `NeuralEnsemble Cookbook <http://neuralensemble.org/cookbook>`_.

If you have questions, there are nice people to answer them at the NeuralEnsemble `googlegroup <http://groups.google.com/group/neuralensemble>`_.

If you feel you have discovered a bug, please let us know so we can fix it by posting a ticket using the "New Ticket" button in the toolbar above.

Usage
-----

.. toctree::
   :maxdepth: 1

   hdf5

API documentation
-----------------

.. toctree::
   :maxdepth: 1

   signals
   spike2
   parameters
   analysis
   stgen
   io
   plotting
   datastore

Developer pages
---------------

For those interested in contributing to the development of NeuroTools, please browse the following links.

.. toctree::
   :maxdepth: 1

   guidelines
   management

Links
=====

-  `http://neuralensemble.org <http://neuralensemble.org/>`_ - The NeuralEnsemble initiative
-  `http://software.incf.org <http://software.incf.org/>`_ - The !INCF software center, a resource that makes it easy for neuroscientists to find, use and share software tools.  -  `NeuralEnsemble Cookbook <http://neuralensemble.org/cookbook>`_ - share your recipies!
-  `PyNN <http://neuralensemble.org/PyNN>`_ - simulator-independent specification of neuronal network models in Python
-  `OpenElectrophy <http://neuralensemble.org/trac/OpenElectrophy>`_ - simplified data and analysis sharing for intra- and extra-cellular recordings leveraging !MySQL.
-  `NEURON <http://www.neuron.yale.edu/neuron>`_ - a standard package in the community for empirically-based simulations of neurons and networks of neurons.
-  `NEST <http://www.nest-initiative.org/>`_ - for simulating large heterogeneous networks of point neurons or neurons with a small number of compartments.
-  `PCSIM <http://sourceforge.net/projects/pcsim/>`_ - a tool for distributed simulation of heterogeneous networks composed of different model neurons and synapses.
-  `FACETS Neuromorphic VLSI <http://wwwasic.kip.uni-heidelberg.de/asicnew/vision/projects/facets/software_and_experiments/>`_ - networks with billions of synapses operating 10\ :sup:`4`\  times faster than biological nervous systems.
-  `Brian <http://brian.di.ens.fr/>`_ - a new simulator for spiking neural networks with in Python, well suited for rapid development of new models, especially networks of single-compartment neurons, and teaching computational neuroscience.
-  `MOOSE/GENESIS <http://moose.sourceforge.net/>`_ - The Multiscale Object-Oriented Simulation Environment. It is the base and numerical core for large, detailed simulations including Computational Neuroscience and Systems Biology.
-  `Neurospaces/GENESIS <http://neurospaces.sourceforge.net/>`_ is a development center for software components of computational neuroscience simulators.
-  `MUSIC <http://incf.org/programs/modeling/music-multi-simulation-coordinator>`_ - MUlti-SImulation Coordinator
-  `PyBrain <http://pybrain.org/>`_ is a modular Machine Learning Library for Python.
-  `BRAHMS <http://brahms.pbwiki.com/>`_ is a Modular Execution Framework, knitting together process software written independently into a single integrated system, and supervising the deployment and execution of that system.
-  `MDP <http://mdp-toolkit.sourceforge.net/>`_ - Modular toolkit for Data Processing (MDP) is a Python data processing framework.
-  `PsychoPy <http://www.psychopy.org/>`_ is an open-source package for creating psychology stimuli in Python combining the graphical strengths of OpenGL with the easy Python syntax.
-  `VisionEgg <http://www.visionegg.org/>`_ is a high level interface between Python and OpenGL for automatic generation of traditional visual stimuli.
-  `PANDORA's Toolbox <http://software.incf.org/software/pandora>`_ is a set of database analysis and visualization tools in MATLAB (c) for simulated and recorded electrophysiological data.


License
=======

NeuroTools:
Analysis, visualization and management of real and simulated neuroscience data. Copyright (C) 2010 Daniel Bruederle, Andrew Davison, Jens Kremkow, Eric Mueller, Eilif Muller, Martin Nawrot, Michael Pereira, Laurent Perrinet, Michael Schmuker, Pierre Yger.

NeuroTools is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

The GNU General Public License does not permit this software to be
redistributed in proprietary programs.

Acknowledgments
===============

NeuroTools is a community-driven open-source project with many hard working contributors: see the ``AUTHORS`` and ``THANKS`` files.

NeuroTools software development is supported in part by the EU under the grant IST-2005-15879 (`FACETS <http://www.facets-project.org/>`_) for years 2005 to 2010 and by  European Union project Number FP7-269921, `BrainScales <http://brainscales.kip.uni-heidelberg.de/>`_ from 2010 to 2014.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

