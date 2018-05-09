# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = Jetstream
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: api test docs

api:
	echo "Building API documents..."
	rm -r api/
	sphinx-apidoc -feMT -o api/ ../jetstream/
	for f in api/*.rst; do\
		perl -pi -e 's/(module|package)$$// if $$. == 1' $$f ;\
	done
	$(MAKE) -C . html

docs: api
	$(MAKE) -C . html

test: api
	sphinx-autobuild . _build

