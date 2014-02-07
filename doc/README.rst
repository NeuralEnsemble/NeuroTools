===========================
Notes on NeuroTools documentation
===========================

The text files in this directory (including this one) are in reStructuredText_ format, which can easily be processed into HTML or LaTeX formats using Docutils_, e.g.::

    $ rst2html.py --initial-header-level=2 parameters.rst > parameters.html

More easily, you may use the provided ```Makefile```:

    $ make html

Many of the files contain examples of interactive python sessions. The validity of this code can be tested by running, for example::

    $ python testdocs.py parameters.rst

.. _reStructuredText: http://docutils.sourceforge.net/rst.html
.. _Docutils: http://docutils.sourceforge.net/