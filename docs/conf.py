"""Sphinx configuration for Thermal Cable Model documentation."""

project = "Thermal Cable Model"
copyright = "2025, Thermal Cable Model Contributors"
author = "Thermal Cable Model Contributors"
release = "0.1.0"
version = "0.1"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_title = "Thermal Cable Model"
html_short_title = "Thermal Cable Model"

autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_class_signature = "separated"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

pygments_style = "friendly"
add_module_names = False
python_use_unqualified_type_names = True
