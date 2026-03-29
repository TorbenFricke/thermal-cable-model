Installation
============

Requirements
------------

Thermal FEM requires **Python 3.9** or later and the following packages:

- `NumPy <https://numpy.org>`_ >= 1.24
- `SciPy <https://scipy.org>`_ >= 1.10
- `Matplotlib <https://matplotlib.org>`_ >= 3.7

Install from source
-------------------

Clone the repository and install in editable (development) mode:

.. code-block:: bash

   git clone https://github.com/your-org/thermal-fem.git
   cd thermal-fem
   pip install -e .

This installs ``thermal_fem`` as an editable package so that changes to the
source code are reflected immediately.

Install dependencies only
-------------------------

If you prefer to run the examples without installing the package, make sure the
dependencies are available:

.. code-block:: bash

   pip install numpy scipy matplotlib

Then add the repository root to your ``PYTHONPATH`` or run scripts from the
repository directory.

Building the documentation
--------------------------

The documentation uses `Sphinx <https://www.sphinx-doc.org>`_ with the
`Read the Docs theme <https://sphinx-rtd-theme.readthedocs.io>`_:

.. code-block:: bash

   pip install sphinx sphinx-rtd-theme
   cd docs
   make html

The built HTML pages are in ``docs/_build/html/``.

Verifying the installation
--------------------------

After installation, verify that the package loads correctly:

.. code-block:: python

   import thermal_fem
   print(thermal_fem.__version__)
   # 0.1.0

You can also run the test suite (requires `pytest <https://pytest.org>`_):

.. code-block:: bash

   pip install pytest
   pytest tests/
