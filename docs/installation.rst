Installation
============

Requirements
------------

Thermal Cable Model requires **Python 3.9** or later and the following packages:

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

This installs ``thermal_cable_model`` as an editable package so that changes to the
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

GitHub Pages (CI deploy)
-------------------------

The repository includes a workflow that builds Sphinx and deploys with
``actions/deploy-pages``. **You must enable Pages once** in the GitHub UI or
deployment returns ``404`` / "Creating Pages deployment failed":

1. Open **Settings → Pages** for the repository.
2. Under **Build and deployment**, set **Source** to **GitHub Actions** (not
   "Deploy from a branch").
3. Re-run the failed workflow (**Actions → Deploy documentation → Re-run jobs**).

**Private repositories:** GitHub Pages for private repos requires a **paid**
plan (e.g. GitHub Pro for personal accounts). On a free account, make the
repository **public** or host docs elsewhere (e.g. Read the Docs).

The site URL is ``https://<owner>.github.io/<repo>/``; ``SPHINX_HTML_BASEURL`` in
the workflow matches this for canonical links.

Verifying the installation
--------------------------

After installation, verify that the package loads correctly:

.. code-block:: python

   import thermal_cable_model
   print(thermal_cable_model.__version__)
   # 0.1.0

You can also run the test suite (requires `pytest <https://pytest.org>`_):

.. code-block:: bash

   pip install pytest
   pytest tests/
