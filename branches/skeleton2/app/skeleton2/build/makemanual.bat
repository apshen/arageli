rem runexamples
latex skeletonman.tex
bibtex skeletonman
latex skeletonman.tex
latex skeletonman.tex
dvips skeletonman.dvi
ps2pdf skeletonman.ps
copy skeletonman.pdf ..\manual.pdf