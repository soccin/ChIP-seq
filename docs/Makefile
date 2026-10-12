#
# Render output.md to output.html and output.pdf. Run on IRIS:
#
#   make -C ChIP-seq/docs
#
# pandoc is ~/bin/pandoc (3.6.4); the PDF engine is the system xelatex
# (TeX Live 2018). lualatex fails there in microtype. Before v1.0.0 the
# HTML came from Marked.app and the PDF from pandoc with LuaLaTeX on a
# laptop, by hand, so the first render from here looks different.
#

PANDOC ?= pandoc

all: output.html output.pdf

output.html: output.md
	$(PANDOC) --standalone --embed-resources \
	    --metadata pagetitle="ChIP-Seq pipeline output" -o $@ $<

output.pdf: output.md
	$(PANDOC) --pdf-engine=xelatex -V geometry:margin=1in -o $@ $<

clean:
	rm -f output.html output.pdf

.PHONY: all clean
