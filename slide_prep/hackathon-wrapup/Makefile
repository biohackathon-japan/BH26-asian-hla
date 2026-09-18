.PHONY: all
all:
	pdflatex -interaction=nonstopmode -halt-on-error slides.tex
	pdflatex -interaction=nonstopmode -halt-on-error slides.tex
	mkdir -p render
	pdftoppm -png -r 140 slides.pdf render/slide
