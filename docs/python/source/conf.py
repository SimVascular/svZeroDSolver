"""Configuration file for the Sphinx documentation builder."""
import svzerodsolver

project = svzerodsolver.NAME
release = svzerodsolver.VERSION
copyright = svzerodsolver.COPYRIGHT

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.napoleon",
    "m2r2",
    "sphinxcontrib.bibtex",
    "sphinx_autodoc_typehints",
]

source_suffix = [".rst", ".md"]
master_doc = "index"
exclude_patterns = []
html_theme = "pydata_sphinx_theme"
# html_logo = "img/logo.png"
# html_favicon = "img/favicon.png"
bibtex_bibfiles = ["refs.bib"]
bibtex_default_style = "unsrt"
set_type_checking_flag = False

html_theme_options = {
    "logo_link": "index",
    "collapse_navigation": True,
    "show_prev_next": False,
    "github_url": "https://github.com/SimVascular/svZeroDSolver",
    "navbar_end": ["navbar-icon-links"],
}
html_title = "%s %s Manual" % (project, release)
