ELSI-examples

This repo contains a few example programs and patches to use ELSI.[1]

The example programs were initiated from ELSI tests files, and modified to read matrices stored in Matrix Market MM format. [2]

The matrices should be symmetric, thus only lower half is considered.

In the MM format considered here each line contains: I, J , VAL(I,J); where I greater than or equal to J. The first line is left for comments.

See Makefile to configure compiler options and ELSI prefix path.

[1] http://elsi-interchange.org/
[2] http://math.nist.gov/MatrixMarket/formats.html

Author: Alvaro Vazquez-Mayagoitia
