Installation
============

Downloading `WDRT`
------------------
`WDRT` is distributed through the `WDRT GitHub repository <https://github.com/WEC-Sim/WDRT/>`_.
The toolbox can either be downloaded via ``git`` or simply by visiting the the `WDRT GitHub repository <https://github.com/WEC-Sim/WDRT/>`_ and downloading it directly.

.. note::

	You may place the downloaded `WDRT` in any location you like; for the remainder of these instructions, the folder that contains the `WDRT` source code will be referred to as ``$WDRT_SOURCE``.

Downloading Python
------------------
If you're new to Python, we suggest one of the following Python distributions:

	* `Anaconda <http://continuum.io/downloads>`_ (Linux, Mac, Windows)

.. note::

	`WDRT` is designed to run with **Python 3** 

Install options
---------------
This page outlines two options for installing the `WDRT`. The first option (`Installing WDRT for Python beginners`_) assumes that you are new to Python and gets you up an running as quickly as possible. The second option (`Installing WDRT for experienced Python users`_) uses a method that allows you to more readily keep you local copy of the `WDRT` update to date with the latest version online.

Installing WDRT for Python beginners
````````````````````````````````````
**Step 1** It can be a pain to install Python, NumPy, SciPy, Matplotlib, h5py and other dependencies that are needed to run `WDRT`. If you're new to Python, the easiest approach is to start by installing the version of either of the following Python distributions:

	* `Anaconda <https://www.anaconda.com/products/individual>`_ (Linux, Mac, Windows)


**Step 2** Download `WDRT` from the `WDRT GitHub repository <https://github.com/WEC-Sim/WDRT/>`_

.. Note:: 

	There may be several branches open on the repository for the purpose of developing new features. 
        The default branch for new users to download is the "master" branch.
        There is no guarantee that the toolbox will work as expected, or at all, if you are working with a non-master branch.

**Step 3** Open a command window (Windows) or terminal window (OSX and Linux) and navigate to ``$WDRT_SOURCE``. 
Once inside the ``$WDRT_SOURCE`` directory execute the following command to install `WDRT`::

	pip install .

**Step 4:** Verify the installation's functionality by running an example located in``$WDRT_SOURCE/examples``

.. code-block:: none

	cd examples
	python example_fatigue.py

Installing WDRT for experienced Python users
````````````````````````````````````````````
The following installation procedure allows for easier updating of the code with ``git``, which can be downloaded `here <https://git-scm.com/downloads>`_.

**Step 1:** Clone or download a copy of `WDRT`::

	git clone https://github.com/WEC-Sim/WDRT $WDRT_SOURCE

**Step 2:** Add the ``$WDRT_SOURCE`` directory to your `PYTHONPATH <https://docs.python.org/3.8/using/cmdline.html?highlight=pythonpath#envvar-PYTHONPATH>`_ `environment variable <https://docs.python.org/3.8/using/cmdline.html?highlight=pythonpath#environment-variables>`_. 

**Step 3:** Verify the installation's functionality by running an example located in``$WDRT_SOURCE/examples``

.. code-block:: none

	cd examples
	python example_fatigue.py

Configuring your PYTHONPATH
````````````````````````````
The following instructions will help you configure your PYTHONPATH, which is a search path Python uses for 
importing other python modules.

**Windows:** 
			| 1.) Navigate to Control Panel -> System -> Advanced system settings -> Environment Variables 
			| 2.) Click "New..." 
			| 3.) Under "Variable Name" type: PYTHONPATH 
			| 4.) Under "Variable Value" enter the location of your Python source directory (i.e "C:\Python")
**Linux/OS X**
			| 1.) Navigate to your home directory
			| 2.) Add a line such as "export PYTHONPATH=“/path/where/your/modules/are/located" to your .bash_rc file if running Linux, or .bash_profile if running OS X
			| 3.) Place the modules you would like to import in the directory you specified in the previous step


Dependencies
-------------
`Python 3 <https://www.python.org/downloads/>`_ and the following Python packages are required to run `WDRT`. 
These packages should be installed using `pip <https://pypi.python.org/pypi/pip>`_  or your preferred package manager:

	* `numpy <http://www.numpy.org>`_
	* `scipy <http://www.scipy.org>`_
	* `matplotlib <http://matplotlib.org>`_
	* `h5py <http://www.h5py.org>`_
	* `sklearn <http://scikit-learn.org/stable/>`_
	* `requests <http://docs.python-requests.org/en/master/>`_
	* `BeautifulSoup4 <https://www.crummy.com/software/BeautifulSoup/>`_
	* `netCDF4  <http://unidata.github.io/netcdf4-python/>`_
	* `statsmodels <http://www.statsmodels.org/>`_
	* `lxml <http://lxml.de/>`_



.. Note::

	the netCDF4 package is only required if you are using a CDIP site in the ESSC module

Troubleshooting
---------------

**Problem:** I can't run any of the examples.

**Solutions:** Check you PYTHONPATH or move the file you want to run into the main WDRT folder.

**Problem:** I can't connect to the NDBC database to download the data I need.

**Solution:** Check your proxy/firewall settings. If you can download data from elsewhere through your proxy/firewall, check the status of the NDBC website with `Down for Everyone <http://downforeveryoneorjustme.com/>`_.

