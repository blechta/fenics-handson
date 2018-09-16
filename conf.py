# -*- coding: utf-8 -*-

project = 'FEniCS hands-on tutorial'
author = 'Jan Blechta, Roland Herzog, Jaroslav Hron, Gerd Wachsmuth'
copyright = '2014, 2015, 2018 ' + author

# The short X.Y version
version = '2017.2'
# The full version, including alpha/beta/rc tags
release = '2017.2.0.dev0'

extensions = [
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinxcontrib.contentui',
    'sphinx.ext.intersphinx',
]

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'README.rst']
pygments_style = 'sphinx'
highlight_language = 'python3'
default_role = 'py:obj'
html_title = 'FEniCS hands-on'
html_theme = 'sphinx_rtd_theme'
html_theme_options = {}
html_static_path = ['_static']

latex_elements = {
    'papersize': 'a4paper',
    'pointsize': '10pt',
    #'preamble': '',
    #'figure_align': 'htbp',
}
latex_documents = [
    (master_doc, 'FEniCS-hands-on.tex', 'FEniCS hands-on', author, 'howto'),
]

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

## Enable ".. only:: solution" sections
#tags.add('solution')

# Fixes for typesetting math
def setup(app):
    app.add_stylesheet('math.css')
    app.add_javascript('copybutton.js')

intersphinx_mapping = {
    'python':('https://docs.python.org/3/', None),
    'ufl': ('https://fenics.readthedocs.io/projects/ufl/en/2017.2.0.post0/', None),
    'dolfin': ('https://fenics.readthedocs.io/projects/dolfin/en/2017.2.0/', None),
    'pydolfinapi': ('https://fenicsproject.org/docs/dolfin/2017.2.0/python/', None),
}
