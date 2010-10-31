#!/bin/sh
../bin/skeleton example.ine
../bin/skeleton cwcv.ine --edges --ineinc
../bin/skeleton cwcv.ext
../bin/skeleton equ.ine --verifyine
../bin/skeleton sf.ine --skeletonformat
../bin/skeleton ucube.ine --avisfukudaformat --adjacency --facetadjacency --extinc --ineinc
../bin/skeleton exvoronoi.ine --ineinc --extinc
../bin/skeleton exdelaunay.ext --extinc
../bin/skeleton ccc6.ext --avisfukudaformat
../bin/skeleton cyc.ine --avisfukudaformat
../bin/skeleton cube.ine
../bin/skeleton cube.ext
../bin/skeleton ex1.ine
../bin/skeleton ex2.ine --skeletonformat
../bin/skeleton ex3.ine --skeletonformat

