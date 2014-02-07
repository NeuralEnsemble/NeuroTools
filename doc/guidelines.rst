Developer Guidelines
====================

Coding Style
------------

-  Please read "The Python Coding Style Guidelines" in
   `PEP-8 <http://www.python.org/dev/peps/pep-0008/>`_ and follow them
   to a reasonable degree.
-  Perhaps use a tool like pylint to help you in improving your coding
   style.

Dependency Checking
-------------------

NeuroTools supplies a dependency checker to gracefully disable functionality
for which the required dependencies are not installed, without affecting other
functionality.  If you would like to use these facilities to add a dependency
to a non-standard python module, please see the "check\_dependency" function or
the various modules that use it for examples.

Testing
-------

If you are looking to contribute functionality, please provide also
tests for said functionality. See the ``test`` folder  for examples of testing.

The complete NeuroTools test suite can be run as follows:

::

    $ nosetests --with-coverage --cover-package=NeuroTools

The python "nose" unittesting suite is a dependency here. On ubuntu/debian,
install it using:

::

    $ sudo apt-get install python-nose python-coverage

Tests for each packages can be run individually, for example for the
``NeuroTools.stgen`` sub-package:

::

    $ python test/test_stgen.py

License of contributed code
---------------------------

All of NeuroTools in under GPL with the copyright remaining with contributors.
In practice, since there maybe many contributions to a given module, we have
until now just attributed a global NeuroTools copyright to all contributors
equally. It is therefore understood that all copyright holders would need vote
to make any licensing changes.  Please contact existing core developers if you
would like to submit code, but have a problem with submitting code to
NeuroTools under these license terms, or if you would like to be added to the
list of NeuroTools copyright holders.

doc strings
-----------

Here are described some guidelines for formatting and content of Python
"doc strings" in NeuroTools

Function doc strings
~~~~~~~~~~~~~~~~~~~~

Here is a template for providing doc strings for functions:

::

    def function(...):
        """
        A description of what is the function doing and returning. Few lines

        Inputs:
        param1 - what is param 1
        param2 - what is param 2
        ....

        Examples:
        >> One or several examples showing how it can be used  

        See also:
        functions related to, used by the function.
        """

Module doc strings
~~~~~~~~~~~~~~~~~~

The top level module doc string should inform the user about the
structure of the module without having to look at the code The module
docstring is the doc you can have if, for example, you do, in ipython
(or with help() in python):

::

    >> import numpy
    >> ?numpy
    >> import scipy.optimize
    >> ?scipy.optimize

As you can see, it's a summary of a module organization, its functions,
and so on... It's a text that should be inserted at the top of the file,
and I suggest the following template:

::

    """
    == Name.of.the.module ==

    A short text describing how the module is useful, crucial, powerful

    Classes
    ---------

    Object1 - The first type of object that can be created with the module
    Object2 - One other type of object and what we can do with


    Functions
    ---------

    function1 - A key function of the module, and what we can do with
    function2 - Same as before...

    """

An example for the stgen module would be as follows:

::

    """
    NeuroTools.stgen
    ================

    A collection of tools for stochastic process generation.

    Classes
    ---------

    StGen - Object to generate stochastic processes of various kinds
            and return them as SpikeTrain or AnalogSignal objects.

    Functions
    -----------

    shotnoise_fromspikes - Convolves the provided spike train with shot
                   decaying exponential.
    gamma_hazard         - Compute the hazard function for a gamma
    process                            with parameters a,b.
    """

