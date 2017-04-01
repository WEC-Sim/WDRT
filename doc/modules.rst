Toolbox Modules
===============
The following sections contain the documentation for `WDRT`.
Note that this documentation is also contained in the python source code.
As such, documentation for any `WDRT` function or classs can also be accessed in iPython via by calling

>>> help(<wdrt class or module>)

WDRT.fatigue module
-------------------------

.. automodule:: WDRT.fatigue
    :members:
    :undoc-members:
    :show-inheritance:

WDRT.shortTermExtreme module
----------------------------------

.. automodule:: WDRT.shortTermExtreme
    :members:
    :undoc-members:
    :show-inheritance:

WDRT.longTermExtreme module
----------------------------------

.. automodule:: WDRT.longTermExtreme
    :members:
    :undoc-members:
    :show-inheritance:


WDRT.ESSC module
----------------------

.. automodule:: WDRT.ESSC
    :members:
    :undoc-members:
    :show-inheritance:

WDRT.MLER module
-------------------------

.. *Currently written in MATLAB and will be converted to Python in the near future*

To improve the model efficiency and reduce required computational time or wave tank testing time, the shipping and offshore industry often perform their mid- and high-fidelity simulations and experimental wave tank tests using an equivalent design wave profile generated from (linear) low-fidelity model solutions. The MLER module was developed based on the most likely extreme response method, which uses the linear RAO of the device and the spectrum for the sea state of interest to produce a focused wave profile that gives the largest response.
This approach does not require an input time series for the device. 
Because it uses linear theory, the underlying assumption is that the higher-order effects are small in comparison to the linear effects.

.. The source code can found in ``$WDRT_SOURCE/WDRT/MLER_toolbox/`` folder and ``testrun.m`` is an example with limited comments.

The source code includes a Python egg in ``$WDRT_SOURCE/WDRT/MLER_toolbox/dist`` for use with ``easy_install``.
