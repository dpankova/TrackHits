Code Review
===========

This project was reviewed by K. Krings, March 30, 2017.

Project information
-------------------

* URL: http://code.icecube.wisc.edu/svn/sandbox/kjero/StartingTrackVeto
* Repository Root: http://code.icecube.wisc.edu/svn
* Repository UUID: 16731396-06f5-0310-8873-f7f720988828
* Revision: 154499
* Last Changed Author: kjero
* Last Changed Rev: 153641
* Last Changed Date: 2017-03-01 04:24:58 +0100 (Wed, 01 Mar 2017)

First glance
------------

The project looks mostly complete: documentation and one example script are
available, unit tests or test scripts are however missing. The example script
expects command line arguments and has no defaults; it should use input file(s)
from our test data if possible.

The code compiles without errors or warnings in both release and debug mode.

Documentation
-------------

The Sphinx documentation is not build when running ``make html`` because there
is no ``DOCS_DIR`` specified in ``CMakeLists.txt``. Moreover, the documentation
contains some syntax errors causing warnings:

    * Line 12: title underline is too short.
    * Lines 15 to 17: blank line and indentation are missing.
    * Lines 20 to 46: same as before plus missing space after \*.
    * Lines 50 to 67: should be re-formatted as a bullet or description list.

The documentation gives a short introduction to what is calculated, describes
the module's output, and which parameters can be set. A more detailed
explanation of the veto is given in the linked presentation. I think, it would
make sense to copy some of the presentation's content to the documentation and
thus make it more detailed. Moreover, I would suggest to remove the reference
to the personal software build, because the path could for example change at
some point or even disappear. Instead, I would change the example in a way that
it is using the specified input files by default.

At least ``StartingTrackVetoUtils.h`` could use some Doxygen documentation,
which should be linked to the Sphinx documentation.

Source code
-----------

Code structure
^^^^^^^^^^^^^^

The project's code structure follows our guidelines.

Coding standards
^^^^^^^^^^^^^^^^

* ``StartingTrackVeto.cxx:36``: parameter has no default value.
* ``StartingTrackVetoUtils.cxx:65``: I would give the OM iterator a more
  meaningful name because the for-loop spans over a lot of lines.
* The header files include libraries that are only used in the implementations
  and should only be included there.
* I have the feeling that the code is using a lot of raw pointers where they
  are actually not needed and could be replaced by references for example.

Readability
^^^^^^^^^^^
In general, the code is easy to understand, but I think the readability would definitely
be improved by introducing more black lines to separate code blocks
that belong together. I also encountered that the usage of tabs and spaces as
well as the code indentation is not always consistent. The code should also be
re-checked for missing or not necessary spaces before and after operators,
respectively; for example ``StartingTrackVetoUtils.cxx:225``.
