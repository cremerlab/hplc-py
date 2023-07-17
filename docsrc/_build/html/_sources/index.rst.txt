
Python Utilities for the Cremer Lab
====================================================
This is a package holding functions useful for parsing, converting, and analyzing 
various types of data produced by members of the Jonas Cremer Lab at Stanford University. 


This is very much in development and will be significantly expanded as the
lab expands! If you have some code that you think others in the group may
want to use (of if you would like some code written for your project) send 
an email to `Griffin Chure <mailto:gchure@stanford.edu>`_.


Installation
------------
If you are interested in contributing to the development, it can be installed 
locally via the `Python Package installer <https://pypi.org/project/pip/>`_ `pip`
via the command line.:: 

   $ pip install cremerlab-utils

If you are interested in contributing to the development of the software, 
you can clone or fork the `GitHub Repository <https://github.com/cremerlab/cremerlab-utils>`_
and install the package locally as follows::

   $ git clone git@github.com:cremerlab/cremerlab-utils
   $ cd cremerlab-utils
   $ pip install -e ./

If you don't want to clone the source repository, you can install the development
version directly from GitHub::

   $ pip install -e git+https://github.com/cremerlab/cremerlab-utils.git

Bug Reports & Questions
-----------------------
The package is licensed under a permissive MIT license and is available on 
`GitHub <https://github.com/cremerlab/cremerlab-utils>`_. If you have questions,
issues, or bugs, please get in touch via the `Git Issues <https://github.com/cremerlab/cremerlab-utils/issues>`_
page.

.. toctree::
   :maxdepth: 3
   :caption: User Guide 

   getting_started/hplc_processing.ipynb
   getting_started/growth_analysis.ipynb

.. toctree::
   :maxdepth: 4
   :caption: API Documentation

   hplc 
   growth 
   bayes


Index
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

