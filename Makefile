# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = FEniCShands-on
SOURCEDIR     = .
BUILDDIR      = _build

# Default build
#all: html
all: html-solution

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: all help Makefile servehtml html-solution upload

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

servehtml:
	browse http://127.0.0.1:8000/_build/html/
	python3 -m http.server --bind 127.0.0.1 2> /dev/null

html-solution:
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" -t solution $(SPHINXOPTS) $(O)

upload: html-solution
	rsync -rP _build/html/ login:public_html/priv/fenics-handson
