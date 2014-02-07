The ``plotting`` module
=======================

This module contains a collection of tools for plotting and image processing that shall
facilitate the generation and handling of your data visualizations.
It utilizes the Matplotlib and the Python Imaging Library (PIL) packages.

Universal Functions and Classes for Normal Matplotlib Use
---------------------------------------------------------

The following functions might be useful for every user of the Matplotlib package.

The function ``get_display(display)``
-------------------------------------

Arguments:
~~~~~~~~~~
* display - if True, a new figure is created. Otherwise, if display is a subplot object, this object is returned.

Returns:
~~~~~~~~
A pylab object with a plot() function to draw the plots.

The function ``progress_bar(progress)``
---------------------------------------

Arguments:
~~~~~~~~~~
* progress - a float between 0. and 1.

Returns:
~~~~~~~~
Prints a progress bar to stdout, filled to the given ratio.

Example of usage::

    >>> progress_bar(0.7)
    |===================================               |

The function ``pylab_params(fig_width_pt, ratio, text_fontsize, tick_labelsize, useTex)``
-----------------------------------------------------------------------------------------

Arguments:
~~~~~~~~~~
* fig_width_pt   - figure width in points. If you want to use your figure inside LaTeX, get this value from LaTeX using '\\showthe\\columnwidth'.
* ratio          - ratio between the height and the width of the figure
* text_fontsize  - size of axes and in-pic text fonts
* tick_labelsize - size of tick label font
* useTex         - enables or disables the use of LaTeX for all labels and texts (for details on how to do that, see http://www.scipy.org/Cookbook/Matplotlib/UsingTex)

Returns:
~~~~~~~~
A dictionary with a set of parameters that help to nicely format figures. 
The return object can be used to update the pylab run command parameters dictionary 'pylab.rcParams'.

The function ``set_axis_limits(subplot, xmin, xmax, ymin, ymax)``
-----------------------------------------------------------------

Arguments:
~~~~~~~~~~
* subplot     - the targeted plot
* xmin, xmax  - the limits of the x axis
* ymin, ymax  - the limits of the y axis

Does:
~~~~~
Defines the axis limits in a plot.

Example of usage::

    >>> x = range(10)
    >>> y = []
    >>> for i in x: y.append(i*i)
    >>> pylab.plot(x,y)
    >>> plotting.set_axis_limits(pylab, 0., 10., 0., 100.)

The function ``set_labels(subplot, xlabel, ylabel)``
----------------------------------------------------

Arguments:
~~~~~~~~~~
* subplot - the targeted plot
* xlabel  - a string for the x label
* ylabel  - a string for the y label

Does:
~~~~~
Defines the axis labels of a plot.

Example of usage::

    >>> x = range(10)
    >>> y = []
    >>> for i in x: y.append(i*i)
    >>> pylab.plot(x,y)
    >>> plotting.set_labels(pylab, 'x', 'y=x^2')

The function ``set_pylab_params(fig_width_pt, ratio, text_fontsize, tick_labelsize, useTex)``
---------------------------------------------------------------------------------------------

Arguments:
~~~~~~~~~~
* fig_width_pt   - figure width in points. If you want to use your figure inside LaTeX, get this value from LaTeX using '\\showthe\\columnwidth'.
* ratio          - ratio between the height and the width of the figure
* text_fontsize  - size of axes and in-pic text fonts
* tick_labelsize - size of tick label font
* useTex         - enables or disables the use of LaTeX for all labels and texts (for details on how to do that, see http://www.scipy.org/Cookbook/Matplotlib/UsingTex)

Does:
~~~~~
Updates a set of parameters within the the pylab run command parameters dictionary 'pylab.rcParams' in order to achieve nicely formatted figures.

Special Plotting Functions and Classes for Specific Requirements
----------------------------------------------------------------

The function ``save_2D_image(mat, filename)``
---------------------------------------------

Arguments:
~~~~~~~~~~
* mat      - a 2D numpy array of floats between 0 and 1
* filename - string specifying the filename where to save the data, has to end on '.png'

Does:
~~~~~
Saves a 2D numpy array of gray shades between 0 and 1 to a PNG file.

Example of usage::

    >>> import numpy
    >>> a = numpy.random.random([100,100]) # creates a 2D numpy array with random values between 0. and 1.
    >>> save_2D_image(a,'randomarray100x100.png')

The function ``save_2D_movie(frame_list, filename, frame_duration)``
--------------------------------------------------------------------

Arguments:
~~~~~~~~~~
* frame_list     - a list of 2D numpy arrays of floats between 0 and 1
* filename       - string specifying the filename where to save the data, has to end on '.zip'
* frame_duration - specifier for the duration per frame, will be stored as additional meta-data for later playing

Does:
~~~~~
Saves a list of 2D numpy arrays of gray shades between 0 and 1 to a zipped tree of PNG files.

Example of usage::

    >>> import numpy
    >>> framelist = []
    >>> for i in range(100): framelist.append(numpy.random.random([100,100])) # creates a list of 2D numpy arrays with random values between 0. and 1.
    >>> save_2D_movie(framelist, 'randommovie100x100x100.zip', 0.1)

The ``SimpleMultiplot`` class
-----------------------------
This class creates a figure consisting of multiple panels, all with the same datatype and the same x-range.

Creation / Constructor Arguments:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* nrows    - number of rows
* ncolumns - number of columns
* title    - the title of the multi-plot
* xlabel   - label for all x-axes
* ylabel   - label for all y-axes
* scaling  - a tuple consisting of two string out of {"liner","log"}, determining the scaling of the x-axis and y-axis

Here is an example of creating a ``SimpleMultiplot`` object::

    >>> nrows = 4
    >>> ncolumns = 5
    >>> title = 'a SimpleMultiplot'
    >>> xlabel = 'the x axis'
    >>> ylabel = 'the y axis'
    >>> scaling = ('linear','log')
    >>> smp = SimpleMultiplot(nrows=self.nrows, ncolumns=self.ncolumns, title=title, xlabel=xlabel, ylabel=ylabel, scaling=scaling)
    
Selecting Panels:
~~~~~~~~~~~~~~~~~
Handles to panels can be directly accessed by their indices via the function call ``panel(i)`` or by stepping iteratively through them with function
``next_panel()``. 

Defining Frames:
~~~~~~~~~~~~~~~~
The frames surrounding a panel can be defined with the function ``set_frame(ax, boollist, linewidth)``, where ax is the handle to the panel of choice,
boollist is a list of four booleans defining if [bottom, left, top, right] of the panel shall have a frame line with width linewidth.

Finalising and Saving:
~~~~~~~~~~~~~~~~~~~~~~
Once a SimpleMultiplot is ready to be saved, calling ``finalise()`` will turn off tick labels for all x-axes except the bottom one.
The whole plot is saved to a filename and type of choice with the call ``save(filename)``.

Autodoc
-------

.. automodule:: NeuroTools.plotting
   :members:
   :undoc-members:

