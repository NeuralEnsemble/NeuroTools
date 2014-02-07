==========
NeuroTools
==========

NeuroTools is a collection of tools for representing and analyzing
neuroscientific data.

NeuroTools is written in Python, and works best with PyNN, or one of the growing
list of simulation engines with a Python front-end such as NEURON, NEST, PCSIM,
FACETS Neuromorphic VLSI, Brian, MOOSE/GENESIS, Neurospaces/GENESIS. NeuroTools
provides modules to facilitate simulation setup, parameterization and
instrumentation, and data management, analysis and visualization. The
data-related tools are equally suited to analysis of experimental data, although
that is not the primary motivation for their development.

NeuroTools aims to:

1. increase the productivity of individual modellers by automating, simplifying,
   and establishing best-practices for common tasks,

2. increase the productivity of the neuroscience modelling community by reducing
   the amount of code duplication across simulation communities,

3. increase the reliability of data analysis tools leveraging Linus's law:
   "given enough eyeballs, all bugs are shallow."

Code status
-----------

.. image:: https://travis-ci.org/NeuralEnsemble/NeuroTools.png?branch=master
   :target: https://travis-ci.org/NeuralEnsemble/NeuroTools
.. image:: https://coveralls.io/repos/NeuralEnsemble/NeuroTools/badge.png
   :target: https://coveralls.io/r/NeuralEnsemble/NeuroTools

This code is considered to be mature and will not evolve - we will  only focus on bug fixes.

Development of a new analysis library takes place @ https://github.com/NeuralEnsemble/elephant

Contributing to NeuroTools
--------------------------

If you are interested in participating in NeuroTools development, please see
the webpage at http://neuralensemble.org/NeuroTools

Related projects
----------------

NeuroTools has evolved in different flavours:

- a legacy version @ https://neuralensemble.kip.uni-heidelberg.de/svn/NeuroTools
- a fork by Trevor Bekolay: https://github.com/tbekolay/neurotools
- another fork (version 0.2.0) : https://github.com/meduz/NeuroTools

Installation
------------

To install use:

    pip install NeuroTools

For installation instructions, see doc/source/install.rst

Documentation
-------------

For documentation see the ``doc`` folder or http://packages.python.org/NeuroTools

For examples see the ``examples`` folder.


More information
----------------

- Home page: http://neuralensemble.org/NeuroTools
- Mailing list: https://groups.google.com/forum/?fromgroups#!forum/neuralensemble
- Documentation: http://packages.python.org/neo/
- Bug reports: https://github.com/NeuralEnsemble/NeuroTools/issues


License
-------

Copyright (C) 2008  Daniel Bruederle, Andrew Davison, Jens Kremkow
Laurent Perrinet, Michael Schmuker, Eilif Muller, Eric Mueller, Pierre Yger

NeuroTools is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

The GNU General Public License does not permit this software to be
redistributed in proprietary programs.

See ``LICENSE``.

