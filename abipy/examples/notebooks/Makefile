# Makefile for the geneation of HTML files.
all: help
html:
	ipython nbconvert --to html *.ipynb
clean: 
	rm *.html

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  clean      to remove standalone HTML files"
