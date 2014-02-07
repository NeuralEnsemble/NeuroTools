Module maintainers
==================

Ideally, each module should have approx. two maintainers.

+-------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| ``signals``:      | Pierre Yger, `LaurentPerrinet <http://invibe.net>`_, Jens Kremkow (others volunteers are more than welcome)                       |
+-------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| ``parameters``:   | Andrew Davison                                                                                                                    |
+-------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| ``io``:           | Pierre Yger and ...                                                                                                               |
+-------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| ``plotting``:     | Daniel Bruederle                                                                                                                  |
+-------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| ``analysis``:     | Eilif Muller                                                                                                                      |
+-------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| ``utilities``:    | Daniel Bruederle                                                                                                                  |
+-------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| ``stgen``:        | Eilif Muller, Michael Schmuker                                                                                                    |
+-------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| ``examples``:     | `LaurentPerrinet <http://invibe.net>`_                                                                                            |
+-------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| ``spike2``:       | Jens Kremkow                                                                                                                      |
+-------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| ``datastore``:    | Andrew Davison                                                                                                                    |
+-------------------+-----------------------------------------------------------------------------------------------------------------------------------+

Maintainers are responsible for

-  identifying missing functionality/tests/docs in their module
-  writing tickets using `GitHub's ticket system <https://github.com/NeuralEnsemble/NeuroTools/issues>`_)
-  finding volunteers to write the code, tests and documentation

**Documentation manager**: Pierre Yger

+-------------------------------+--------------------+
| advertising:                  | Laurent Perrinet   |
+-------------------------------+--------------------+

*Responsibilities*: combining the documentation from the different modules into
a coherent whole, ensuring consistent formatting, spell-checking, etc.

**Testing manager**: Andrew Davison

*Responsibilities*:

-  identifying areas of the codebase that are not well tested, and notifying
  the module maintainers
-  organizing/collecting tests that use several of the NeuroTools modules, i.e.
   integration tests rather than unit tests.

**Packaging manager**: Eric Mueller

*Responsibilities*:

-  setup.py, i.e. making sure that distutils installation works.
-  uploading packages to PyPI, software.incf.org, etc.
-  evaluating whether easy\_install would work for NeuroTools
