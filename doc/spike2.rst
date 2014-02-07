The ``spike2`` module
=========================

The ``spike2`` module offers an easy way for reading data from CED's Spike2 Son files into the NeuroTools environment.

The main interaction with with Spike2 files is provided by the SON Library which was written by::

 Antonio Gonzalez
 Department of Neuroscience
 Karolinska Institutet
 Antonio.Gonzalez at cantab.net
 http://www.neuro.ki.se/broberger/

The ``spike2`` module integrates the SON library such that the loading function directly returns ``NeuroTools.signals`` objects.
This makes it very easy to apply certain analysis routines on simulated and real experimental data.

Loading channels
----------------

Usually SON files contain data from multiple channels. The channels can contain analog signals, like membrane potential traces or discrete data like markers. The ``spike2`` module knows which type of channel you are loading and thus returns the appropriate ``NeuroTools.signals`` object.

Here some examples how one can simply load the data. For the following example we assume that the SON file contains the membrane potential trace in channel 1, the markers in channel 2 and in channel 3 is data that we dont want to use:

To load the ``spike2`` module you could do the following::

	>>> import NeuroTools.spike2.spike2channels as spike2

to load all channels in the file just use the load function::

	>>> all_channels = spike2.load(filename)
	>>> all_channels.keys()
	[1,2,3]

since we don't need channel 3 we can use the channels parameter to specify which channel to load::

	>>> just_needed_channels = spike2.load(filename, channels=[1,2])
	>>> just_needed_channels.keys()
	[1,2]

what objects are returned?::

	>>> type(just_needed_channels[1])
 	<class 'NeuroTools.spike2.spike2channels.Adc'>
	>>> type(just_needed_channels[2])
	<class 'NeuroTools.spike2.spike2channels.Marker'>

what can I now do with it? I.e. the Adc channel type is inherited from ``signals.AnalogSignal`` and thus you could plot an event_triggered_average on the marker times (please note that the markers are in seconds, we need them in milliseconds therefore the factor 1000)::

	>>> vm = just_needed_channels[1]
	>>> marker = just_needed_channels[2]
	>>> vm..event_triggered_average(marker.times*1000., display=True)

in case your vm channel contains multiple stimuli conditions that you dont want to average, but still be separated you can use the slice_by_events function::

	>>> vm_sclices = vm.slice_by_events(marker.times*1000.,t_min=100,t_max=1000)

For further examples of what to do with ``NeuroTools.signals`` objects, please refer to the documentation for the ``signals`` module.

Autodoc
-------

.. automodule:: NeuroTools.spike2
   :members:
   :undoc-members:


