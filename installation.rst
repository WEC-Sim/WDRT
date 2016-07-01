Installation
==========================

Downloading `WDRT`
------------------------
`WDRT` is distributed through the `WDRT GitHub repository <https://github.com/WEC-Sim/WDRT/>`_. The toolbox can either be downloaded via ``git`` or simply by visiting the the `WDRT GitHub repository <https://github.com/WEC-Sim/WDRT/>`_ and downloading it directly.

.. note::

	You may place the downloaded `WDRT` in any location you like; for the remainder of these instructions, the folder that contains the `WDRT` source code will be referred to as ``$WDRT_SOURCE``.

Dependencies
-------------
`Python 2.7.x <https://www.python.org/downloads/>`_ and the following Python packages are required to run `WDRT`.
These packages can easily be installed using using `pip <https://pypi.python.org/pypi/pip>`_  or your preferred package installation method:

	* `numpy <http://www.numpy.org>`_
	* `scipy <http://www.scipy.org>`_
	* `matplotlib <http://matplotlib.org>`_
	* `h5py <http://www.h5py.org>`_
	* `sklearn <http://scikit-learn.org/stable/>`_
	* `requests <http://docs.python-requests.org/en/master/>`_
	* `BeautifulSoup4 <https://www.crummy.com/software/BeautifulSoup/>`_

Alternatively, if you're new to Python, these packages can be obtained using one of following Python distributions:

	* `Python(x,y) <http://python-xy.github.io>`_
	* `Anaconda <https://www.continuum.io/downloads>`_
	* `Canopy <https://www.enthought.com/products/canopy/>`_

Installation Process
--------------------
The following installation procedure allows for easier updating of the code with ``git``

**Step 1:** Clone or download a copy of `WDRT`::

	git clone http://wec-sim.github.io/WDRT $WDRT_SOURCE

**Step 2:** Add the ``$WDRT_SOURCE`` directory to your `PYTHONPATH <https://docs.python.org/2/using/cmdline.html#environment-variables>`_ environment variable (`Windows <https://docs.python.org/2/using/windows.html#excursus-setting-environment-variables>`_, `Mac OSX <https://docs.python.org/2/using/mac.html?highlight=pythonpath#configuration>`_, `Linux <https://wiki.archlinux.org/index.php/Environment_variables>`_).

**Step 3:** Verify the installation's functionality by running the examples located in ``$WDRT_SOURCE/examples``.
