.. image:: https://upload.wikimedia.org/wikipedia/commons/thumb/1/1c/USGS_logo_green.svg/320px-USGS_logo_green.svg.png
        :target: http://www.usgs.gov/
        :alt: U.S. Geological Survey logo

waterpy
===============================

*waterpy* is a rainfall-runoff model that predicts the amount of water
flow in rivers. *waterpy* is a command line application written in Python
using Click_, and is a complete conversion of the original rainfall-runoff
model, called Topmodel, from Fortran to Python. The specific version of 
Topmodel that *waterpy* was initially based on a version by David Wolock, 
U.S. Geological Survey. Please see report below for more details: 

        Wolock, D.M., "Simulating the variable-source-area concept of
        streamflow generation with the watershed model Topmodel", U.S. Geological
        Survey, Water-Resources Investigations Report 93-4124, 1993.

The David Wolock version was recoded and reproduced in the `topmodelpy`_ project.

*waterpy* is a fork of the David Wolock version with many modifications in an 
attempt to replicate the Topmodel versions by Leon Kaufmann (USGS) and
Tanja Williamson (USGS). Please see report below for more details:

        Williamson, T.N., Lant, J.G., Claggett, P.R., Nystrom, E.A.,
        Milly, P.C.D., Nelson, H.L., Hoffman, S.A., Colarullo, S.J., and Fischer, J.M.,
        2015, Summary of hydrologic modeling for the Delaware River Basin using the
        Water Availability Tool for Environmental Resources (WATER): U.S. Geological
        Survey Scientific Investigations Report 2015–5143, 68 p.,
        http://dx.doi.org/10.3133/sir20155143


Features
--------

* Written entirely in Python for ease of use and model extension.
* Generates a report.hmtl file of model results with interactive plots.


Example
-------

<<<<<<< HEAD
To run the provided example edit the provided .ini files and insert the path names for the
input/output directory and location of the database directory.
=======
>>>>>>> b21aa33c4be91bdba14a980db6e5ac9e0ce56529
To run waterpy, give waterpy the command `run` along with the path to the 
model configuration file:

::

    $ waterpy run <path-to-your-modelconfig.ini>
    $ waterpy run data/modelconfig.ini

The model results are saved to the output directory location specified in
the model configuration file.

A sample model configuration file called `modelconfig.ini` is located in the 
`data/` directory along with sample input files located in the `inputs/`
directory and sample output files located in the `outputs/` directory.

<<<<<<< HEAD
The use of the geospatial tools requires a download/creationg of relevant files.
A database compiled for the state of Kentucky is available for download from
ScienceBase at https://www.sciencebase.gov/catalog/item/5f7e339682ce1d74e7dda351.
Geospatial tools also require a polygon shapefile of the basin area.  An example
shapefile is provided for use with this database in /example/Hope
The geospatial tools can be run from the command line.  Nagivate to the sample
directory as geospatial.py:::

    $ python geospatial.py

Requirements
------------

This packages requires python3.  Use of the geospatial tools requires Python 3.7, 3.8 or 3.9
if packages are installed through pip install.
The following are the main requirements/dependencies:

astroid==2.11.2

certifi==2020.11.8

cftime==1.3.0

click==7.1.2

colorama==0.4.4

configparser==5.0.1

cycler==0.10.0

dill==0.3.4

GDAL==3.1.2

isort==5.10.1

Jinja2==2.11.2

kiwisolver==1.3.1

lazy-object-proxy==1.7.1

MarkupSafe==1.1.1

matplotlib==3.3.3

mccabe==0.7.0

mpld3==0.5.1

netCDF4==1.5.4

numpy==1.19.4

pandas==1.1.4

pathlib==1.0.1

Pillow==8.0.1

platformdirs==2.5.2

PyCRS==1.0.2

pyparsing==2.4.7

pyproj==3.0.0.post1

python-dateutil==2.8.1

pytz==2020.4

scipy==1.5.4

Shapely==1.7.1

six==1.15.0

tomli==2.0.1

typed-ast==1.5.3

typing-extensions==4.2.0

wrapt==1.14.0 
=======

Documentation
-------------

An initial website of documentation is located in the `docs` directory.


Tests
-----

A suite of tests were built using `pytest <http://pytest.org/latest/>`_.

To run the test suite, from the command line in the project's root directory::

    $ py.test tests/



Requirements
------------

The following are the main requirements/dependencies:

Click==7.0    

Jinja2==2.10.1      

matplotlib==3.0.3     

mpld3==0.3     

numpy==1.16.2     

pandas==0.24.1     

pytest==4.3.0     

scipy==1.2.1     
>>>>>>> b21aa33c4be91bdba14a980db6e5ac9e0ce56529


Installation
------------

To install waterpy from source:

1. Check that you have Python_ installed::

    $ python --version

If you do not have Python_ installed, please download the latest version from `Python's download page`_

2. Download waterpy from the repository and extract to a directory of your choice.

   Or, if you have git_ installed you can clone the project::

    $ git clone <remote url to waterpy>

3. Navigate to the project's root directory where the setup script called `setup.py` is located::

    $ cd waterpy/

| The `setup.py` is a Python file that contains information regarding the installation of a Python module/package, and
| usually specifies that the module/package has been packaged and distributed with the standard Python distribution
| package called Distutils_.

<<<<<<< HEAD
4. To install GDAL python geospatial .whl files are included for Python 3.7-3.9.::  

	$ pip install /utilities/GDAL-3.1.2-cp38-cp38-win_amd64.whl

	replace cp38 with you version of python.

5. pip install requirements.txt
	
6. Run `setup.py` with the `install` command::

    $ python setup.py install 

waterpy will now be installed to the standard location for third-party Python modules on your
computer platform.
=======
4. Run `setup.py` with the `install` command::

    $ python setup.py install

waterpy will now be installed to the standard location for third-party Python modules on your
computer platform.

All the required 3rd party packages in the "Requirements" section should be installed as well,
however, if the additional packages are not installed you can manually install them using `pip`::

    $ pip install click jinja2 matplotlib mpld3 numpy pandas scipy

or from within the parent *waterpy* directory::

    $ pip install requirements.txt
>>>>>>> b21aa33c4be91bdba14a980db6e5ac9e0ce56529

For more information regarding installing third-party Python modules, please see `Installing Python Modules`_
For a description of how installation works including where the module will be installed on your computer platform,
please see `How Installation Works`_.


License
-------

This software is licensed under `CC0 1.0`_ and is in the `public domain`_ because it contains materials that originally
came from the `U.S. Geological Survey (USGS)`_, an agency of the `United States Department of Interior`_. For more
information, see the `official USGS copyright policy`_.

.. image:: http://i.creativecommons.org/p/zero/1.0/88x31.png
        :target: http://creativecommons.org/publicdomain/zero/1.0/
        :alt: Creative Commons logo


Disclaimer
----------

<<<<<<< HEAD
This software has been approved for release by the U.S. Geological Survey (USGS). Although the software has been
subjected to rigorous review, the USGS reserves the right to update the software as needed pursuant to further analysis
and review. No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of 
the software and related material nor shall the fact of release constitute any such warranty. Furthermore, the 
software is released on condition that neither the USGS nor the U.S. Government shall be held liable for any damages 
resulting from its authorized or unauthorized use.
=======
This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely
best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed
or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor
shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS
nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the
software.
>>>>>>> b21aa33c4be91bdba14a980db6e5ac9e0ce56529

The USGS provides no warranty, expressed or implied, as to the correctness of the furnished software or the suitability
for any purpose. The software has been tested, but as with any complex software, there could be undetected errors. Users
who find errors are requested to report them to the USGS.

References to non-USGS products, trade names, and (or) services are provided for information purposes only and do not
constitute endorsement or warranty, express or implied, by the USGS, U.S. Department of Interior, or U.S. Government, as
to their suitability, content, usefulness, functioning, completeness, or accuracy.

Although this program has been used by the USGS, no warranty, expressed or implied, is made by the USGS or the United
States Government as to the accuracy and functioning of the program and related program material nor shall the fact of
distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.

This software is provided "AS IS."


Author
------

Alex Headman <AHeadman@usgs.gov>
Jeremiah Lant <jlant@usgs.gov>


.. _Python: https://www.python.org/
<<<<<<< HEAD
.. _U.S. Geological Survey: https://www.usgs.gov/
.. _United States Department of Interior: https://www.doi.gov/
.. _official USGS copyright policy: http://www.usgs.gov/visual-id/credit_usgs.html#copyright/
.. _U.S. Geological Survey (USGS) Software User Rights Notice: http://water.usgs.gov/software/help/notice/
=======
.. _pytest: http://pytest.org/latest/
.. _Click: https://click.palletsprojects.com/
.. _Sphinx: http://sphinx-doc.org/
.. _public domain: https://en.wikipedia.org/wiki/Public_domain
.. _CC0 1.0: http://creativecommons.org/publicdomain/zero/1.0/
.. _U.S. Geological Survey: https://www.usgs.gov/
.. _USGS: https://www.usgs.gov/
.. _U.S. Geological Survey (USGS): https://www.usgs.gov/
.. _United States Department of Interior: https://www.doi.gov/
.. _official USGS copyright policy: http://www.usgs.gov/visual-id/credit_usgs.html#copyright/
.. _U.S. Geological Survey (USGS) Software User Rights Notice: http://water.usgs.gov/software/help/notice/
.. _Python's download page: https://www.python.org/downloads/
.. _git: https://git-scm.com/
.. _Distutils: https://docs.python.org/3/library/distutils.html
.. _Installing Python Modules: https://docs.python.org/3.5/install/
.. _How Installation Works: https://docs.python.org/3.5/install/#how-installation-works
.. _topmodelpy: https://github.com/jlant/topmodelpy
>>>>>>> b21aa33c4be91bdba14a980db6e5ac9e0ce56529
