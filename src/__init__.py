# -*- coding: utf-8 -*-
"""
NeuroTools
==========

NeuroTools is a collection of tools for
representing and anlyzing neuroscientific data.

For more information see:
http://neuralensemble.org/NeuroTools

Available subpackages
---------------------

signals
    Provides core classes for manipulation of spike trains and analog signals.
spike2
    Offers an easy way for reading data from CED's Spike2 Son files.
parameters
    Contains classes for managing large, hierarchical parameter sets.
analysis
    Cross-correlation, tuning curves, frequency spectrum, etc.
stgen
    Various stochastic process generators relevant for Neuroscience
    (OU, poisson, inhomogenous gamma, ...).
io
    NeuroTools support for reading and writing of files in various formats.
plotting
    Routines for plotting and visualization.
datastore
    A consistent interface for persistent data storage
    (e.g. for caching intermediate results).

Subpackage specific documentation is available by importing the
subpackage, and requesting help on it::

  >>> import NeuroTools.signals
  >>> help(NeuroTools.signals)
  ... # doctest: +SKIP

Utilities
---------
NeuroTools contains some fancy logging and dependency checking mechamisms,
at least for now (might remove them soon).

"""

__all__ = ['analysis', 'parameters', 'plotting', 'signals', 'stgen',
           'io', 'datastore', 'spike2', 'tisean']
__version__ = "0.2.0 (Bursting Basketcell)"
import warnings


#########################################################
## ALL DEPENDENCIES SHOULD BE GATHERED HERE FOR CLARITY
#########################################################

# The nice thing would be to gather every non standard
# dependency here, in order to centralize the warning
# messages and the check
dependencies = {'pylab': {'website': 'http://matplotlib.sourceforge.net/',
                           'is_present': False, 'check': False},
                'matplotlib': {'website': 'http://matplotlib.sourceforge.net/',
                               'is_present': False, 'check': False},
                'tables': {'website': 'http://www.pytables.org/moin',
                           'is_present': False, 'check': False},
                'PIL': {'website': 'http://www.pythonware.com/products/pil/',
                        'is_present': False, 'check': False},
                'scipy': {'website': 'http://numpy.scipy.org/',
                          'is_present': False, 'check': False},
                'rpy': {'website': 'http://rpy.sourceforge.net/',
                        'is_present': False, 'check': False},
                'rpy2': {'website': 'http://rpy.sourceforge.net/rpy2.html',
                         'is_present': False, 'check': False},
                'IPython': {'website': 'http://ipython.scipy.org/',
                            'is_present': False, 'check': False},
                'interval': {'website': 'http://pypi.python.org/pypi/interval/1.0.0',
                            'is_present': False, 'check': False},
                'TableIO' : {'website': 'http://kochanski.org/gpk/misc/TableIO.html',
                                'is_present': False, 'check': False},
                ## Add here your extensions ###
               }

#########################################################
## Function to display error messages on the dependencies
#########################################################


class DependencyWarning(UserWarning):
    pass

def get_import_warning(name):
    return '''** {} ** package is not installed.
To have functions using {} please install the package.
website : {}
'''.format(name, name, dependencies[name]['website'])

def get_runtime_warning(name, errmsg):
    return """** {} ** package is installed but cannot be imported.
              The error message is: {}""".format(name, errmsg)

def check_numpy_version():
    import numpy
    numpy_version = numpy.__version__.split(".")[0:2]
    numpy_version = float(".".join(numpy_version))
    if numpy_version >= 1.2:
        return True
    else:
        return False

def check_pytables_version():
    import tables
    if tables.__version__ <= 2:
        raise Exception("""PyTables version must be >= 1.4,
                        installed version is {}""".format(__version__))

def check_dependency(name):
    if dependencies[name]['check']:
        return dependencies[name]['is_present']
    else:
        try:
            exec("import {}".format(name))
            dependencies[name]['is_present'] = True
        except ImportError:
            warnings.warn(get_import_warning(name), DependencyWarning)
        except RuntimeError as errmsg:
            warnings.warn(get_runtime_warning(name, errmsg), DependencyWarning)
        dependencies[name]['check'] = True
        return dependencies[name]['is_present']

# Setup fancy logging
# red = 0010
# green = 0020
# yellow = 0030
# blue = 0040
# magenta = 0050
# cyan = 0060
# bright = 0100

try:
    import ll.ansistyle

    def colour(col, text):
        try:
            return unicode(ll.ansistyle.Text(col, unicode(text)))
        except UnicodeDecodeError as e:
            raise UnicodeDecodeError("{}. text was {}".format(e, text))
except ImportError:

    def colour(col, text):
        return text

import logging

# Add a header() level to logging (higher than warning, lower than error)
logging.HEADER = (logging.WARNING + logging.ERROR) / 2
logging.addLevelName(logging.HEADER, 'HEADER')

root = logging.getLogger()


def root_header(msg, *args, **kwargs):
    if len(root.handlers) == 0:
        logging.basicConfig()
    apply(root.header, (msg,) + args, kwargs)


def logger_header(self, msg, *args, **kwargs):
    if self.manager.disable >= logging.HEADER:
        return
    if logging.HEADER >= self.getEffectiveLevel():
        apply(self._log, (logging.HEADER, msg, args), kwargs)

logging.Logger.header = logger_header
logging.header = root_header


class FancyFormatter(logging.Formatter):
    """A log formatter that colours and indents
    the log message depending on the level.

    """

#     DEFAULT_COLOURS = {
#         'CRITICAL': bright + red,
#         'ERROR': red,
#         'WARNING': magenta,
#         'HEADER': bright + yellow,
#         'INFO': cyan,
#         'DEBUG': green
#     }

    DEFAULT_INDENTS = {
        'CRITICAL': "",
        'ERROR': "",
        'WARNING': "",
        'HEADER': "",
        'INFO': "  ",
        'DEBUG': "    ",
    }

    def __init__(self, fmt=None, datefmt=None,
#                  colours=DEFAULT_COLOURS,
                 mpi_rank=None):
        logging.Formatter.__init__(self, fmt, datefmt)
        self._colours = colours
        self._indents = FancyFormatter.DEFAULT_INDENTS
        if mpi_rank is None:
            self.prefix = ""
        else:
            self.prefix = str(mpi_rank)

    def format(self, record):
        s = logging.Formatter.format(self, record)
        if record.levelname == "HEADER":
            s = "=== {} ===".format(s)
        if self._colours:
            s = colour(self._colours[record.levelname], s)
        return self.prefix + self._indents[record.levelname] + s


class NameOrLevelFilter(logging.Filter):
    """Logging filter which allows messages that either have an approved name,
    or have a level >= the level specified.

    The intended use is when you want to receive most messages at a high level,
    but receive certain named messages at a lower level, e.g. for debugging a
    particular component.

    """
    def __init__(self, names=[], level=logging.INFO):
        self.names = names
        self.level = level

    def filter(self, record):
        if len(self.names) == 0:
            allow_by_name = True
        else:
            allow_by_name = record.name in self.names
        allow_by_level = record.levelno >= self.level
        return (allow_by_name or allow_by_level)


def init_logging(filename, file_level=logging.INFO,
                 console_level=logging.WARNING, mpi_rank=None):
    if mpi_rank is None:
        mpi_fmt = ""
    else:
        mpi_fmt = str(mpi_rank)
    logging.basicConfig(level=file_level,
                        format='%%(asctime)s {}%%(name)-10s %%(levelname)-6s' +
                          '%%(message)s [%%(pathname)s:%%(lineno)d]'.format(mpi_fmt),
                        filename=filename,
                        filemode='w')
    console = logging.StreamHandler()
    console.setLevel(console_level)
    console.setFormatter(FancyFormatter('%(message)s', mpi_rank=mpi_rank))
    logging.getLogger('').addHandler(console)
    return console