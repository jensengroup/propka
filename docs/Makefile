# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile ghpages

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# push to this repo's gh-pages (must be named origin and must have a gh-pages branch)

REPONAME = origin

ghpages: html
	-(ORIGIN_REMOTE_URL=`git config --get remote.$(REPONAME).url`; \
	rev=`git rev-parse --short HEAD`; \
        cd $(BUILDDIR)/$^ \
	&& git init \
	&& git remote add origin $${ORIGIN_REMOTE_URL} \
	&& git fetch --depth 50 origin gh-pages \
	&& git reset origin/gh-pages \
	&& touch .nojekyll \
	&& git add -A . \
        && git diff-index --quiet HEAD -- || git commit -m "updated sphinx docs from $${rev}" \
	&& git push -q origin HEAD:gh-pages)
	rm -rf $(BUILDDIR)/$^/.git
