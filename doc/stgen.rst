The ``stgen`` module
====================

This module offers various stochastic generators for point processes that can
be used as spike trains. 

The StGen class
---------------

Creation
~~~~~~~~

Create an ``StGen`` object:

    >>> st_gen = StGen()

This will initialize the stochastic generator and by default try to create a
numpy random generator instance. 
    
Optionally, you can also pass a random number generator instance to the 
constructor:

    >>> import numpy
    >>> st_gen = StGen(rng = numpy.random.RandomState())

You can also use random number generators from gnu scientific library (gsl):

    >>> from pygsl.rng import rng
    >>> st_gen_gsl = StGen(rng = rng())

If you want to seed the random number generator with a specific seed, you can 
do so in the constructor:

    >>> st_gen = StGen(seed = 1234567)

Alternatively, you can re-seed the random number generator when the StGen 
object has already been created:

    >>> st_gen.seed(7654321)
    
Poisson-distributed point processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the ``StGen``-object, you can generate point processes with 
inter-spike-intervals distributed according to a poisson distribution:

    >>> st_gen = StGen()
    >>> spike_train_poisson = st_gen.poisson_generator(rate = 100., 
                                                       tstart = 0., 
                                                       tstop = 2500.)

This generates a NeuroTools.SpikeTrain object, containing spike times with an 
approximate rate of 100 Hz and a duration of 2.5 seconds. 

If you want a numpy array of spike times rather than a SpikeTrain object, 
specify the array keyword:

    >>> spike_train_array = st_gen.poisson_generator(rate = 100., array = True)

Dynamic poisson-distributes point processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

StGen can also generate inhomogeneous poisson processes, i.e. spike trains with 
dynamically changing rates:
    
    >>> spike_train_dyn = st_gen.poissondyn_generator(rate = [50., 80., 30.], 
                                                      t = [0., 1000., 2000.], 
                                                      tstop = 2.5,
                                                      array = False)

This will generate a SpikeTrain object containing spike times with an 
approximate rate of 50 Hz for one second, followed by 80 Hz for one second, and 
finally 30 Hz for half a second. Note that t[0] is used as tstart.

Autodoc
-------

.. automodule:: NeuroTools.stgen
   :members:
   :undoc-members:


