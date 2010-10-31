..\bin\skeleton32f.exe example.ine
..\bin\skeleton32f.exe cwcv.ine --edges --ineinc
..\bin\skeleton32f.exe cwcv.ext
..\bin\skeleton32f.exe equ.ine --verifyine
..\bin\skeleton32f.exe sf.ine --skeletonformat
..\bin\skeleton32f.exe ucube.ine --avisfukudaformat --adjacency --facetadjacency --extinc --ineinc
..\bin\skeleton32f.exe exvoronoi.ine --ineinc --extinc
..\bin\skeleton32f.exe exdelaunay.ext --extinc

..\bin\skeleton32f.exe ccc6.ext --avisfukudaformat
..\bin\skeleton32f.exe cyc.ine --avisfukudaformat

..\bin\skeleton32f.exe cube.ine
..\bin\skeleton32f.exe cube.ext

..\bin\skeleton32f.exe ex1.ine
..\bin\skeleton32f.exe ex2.ine --skeletonformat
..\bin\skeleton32f.exe ex3.ine --skeletonformat

fc golden.out\* *



