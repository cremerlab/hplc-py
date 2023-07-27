# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'hplc-py'
copyright = '2023, Griffin Chure & Jonas Cremer'
author = 'Griffin Chure & Jonas Cremer'
release = '0.0.01'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.napoleon', 'nbsphinx', 'sphinx.ext.coverage', 
              'sphinx.ext.autodoc']

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = "_static/logo_horizontal-02.png"
html_theme_options = {
    "logo_only": True,
    "sticky_navigation": True,
    "collapse_navigation": True,
    "style_nav_header_background":"#3C3E47",
}

html_css_files = ['css/custom.css']
# html_theme_options = {
#     "logo": {
#         "image_dark": '_static/logo_horizontal-02.svg',
#         "image_light": '_static/logo_horizontal-01.svg',
#     } ,
#     "navbar_end": ["navbar-icon-links"]
#     }   
# html_context = {"default_mode": "light"}


# html_css_files = [
#     'css/custom.css',
# ]