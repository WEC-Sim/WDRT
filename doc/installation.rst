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
	* `PythonXY <https://code.google.com/p/pythonxy/>`_ (Windows)

.. note::

	`WDRT` is currently designed to run with **Python 2.7.x**. `WDRT` currently will not fully run on Python 3.5.x. 

Install options
---------------
This page outlines two options for installing the `WDRT`. The first option (`Installing WDRT for Python beginners`_) assumes that you are new to Python and gets you up an running as quickly as possible. The second option (`Installing WDRT for experienced Python users`_) uses a method that allows you to more readily keep you local copy of the `WDRT` update to date with the latest version online.

Installing WDRT for Python beginners
````````````````````````````````````
**Step 1** It can be a pain to install Python, NumPy, SciPy, Matplotlib, h5py and other dependencies that are needed to run `WDRT`. If you're new to Python, the easiest approach is to start by installing the 2.7.x version of either of the following Python distributions:

	* `Anaconda <http://continuum.io/downloads>`_ (Linux, Mac, Windows)
	* `PythonXY <https://code.google.com/p/pythonxy/>`_ (Windows)

.. Note::

	If you are using Anaconda, PythonXY or another non-Python.org distribution, you may need to install one of the packages identified in the `Dependencies`_ section. If you need to do so, you should install any needed modules via your distribution's package manager.

**Step 2** Download `WDRT` from the `WDRT GitHub repository <https://github.com/WEC-Sim/WDRT/>`_
.. Note:: 

	There may be several branches open on the repository at any given time for the purpose of developing new features. ONLY the Master branch is meant to be used outside of our development team. There is no guarantee that the toolbox will work as expected, or at all, if you are working with a non-master branch.

**Step 3** Open a command window (Windows) or terminal window (OSX and Linux) and navigate to ``$WDRT_SOURCE``. Once inside the ``$WDRT_SOURCE`` directory execute the following command to install `WDRT`::

	python setup.py install --user

**Step 4:** Verify the installation's functionality by running the examples located in``$WDRT_SOURCE/examples``

.. code-block:: none

	cd examples
	python example_envSamplying.py
	python example_contourApproach.py
	python example_shortTermExtreme.py
	python example_fatigue.py

Installing WDRT for experienced Python users
````````````````````````````````````````````
The following installation procedure allows for easier updating of the code with ``git``, which can be downloaded `here <https://git-scm.com/downloads>`_.

**Step 1:** Clone or download a copy of `WDRT`::

	git clone http://wec-sim.github.io/WDRT $WDRT_SOURCE

**Step 2:** Add the ``$WDRT_SOURCE`` directory to your `PYTHONPATH <https://docs.python.org/2/using/cmdline.html#environment-variables>`_ environment variable (`Windows <https://docs.python.org/2/using/windows.html#excursus-setting-environment-variables>`_, `Mac OSX <https://docs.python.org/2/using/mac.html?highlight=pythonpath#configuration>`_, `Linux <https://wiki.archlinux.org/index.php/Environment_variables>`_). 

**Step 3:** Verify the installation's functionality by running the examples located in``$WDRT_SOURCE/examples``

.. code-block:: none

	cd examples
	python example_envSamplying.py
	python example_contourApproach.py
	python example_shortTermExtreme.py
	python example_fatigue.py

Dependencies
-------------
`Python 2.7.x <https://www.python.org/downloads/>`_ and the following Python packages are required to run `WDRT`. `WDRT` currently will not fully run on Python 3.5.x.
These packages can easily be installed using using `pip <https://pypi.python.org/pypi/pip>`_  or your preferred package installation method:

	* `numpy <http://www.numpy.org>`_
	* `scipy <http://www.scipy.org>`_
	* `matplotlib <http://matplotlib.org>`_
	* `h5py <http://www.h5py.org>`_
	* `sklearn <http://scikit-learn.org/stable/>`_
	* `requests <http://docs.python-requests.org/en/master/>`_
	* `BeautifulSoup4 <https://www.crummy.com/software/BeautifulSoup/>`_

Troubleshooting
---------------

**Problem:** I can't run any of the examples.

**Solutions:** Check you PYTHONPATH or move the file you want to run into the main WDRT folder.

**Problem:** I can't connect to the NDBC database to download the data I need.

**Solution:** Check your proxy/firewall settings. If you can download data from elsewhere through your proxy/firewall, check the status of the NDBC website with `Down for Everyone <http://downforeveryoneorjustme.com/>`_.

**Problem:** I want to use the MLER toolbox, but it's in an .egg file. 

**Solution:** Make sure you have the easy install package, `which can be downloaded here <https://pypi.python.org/pypi/setuptools>`_. Then, run the following command in the command line::

	Python -m easy_install C:\path\to\mler.egg
