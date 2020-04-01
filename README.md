MDreac
======

Molecular Dynamics simulation coupled with chemical reactions.

In the directory `oldsrc` you find my original MDreac software package. I developed the package
during my graduate studies, and the algorithms and results are described in my Ph.D. thesis. It
is written in Fortran-77 (with a few common extensions), and it is possible to compile it using
GNU Fortran v8 or later.

A parallel version is included. It is done using PVM.

md2.c
-----
A newer - and simpler - version is located in the directory `src`. It is written in C and was
developed as part of a series of feature articles on multicore programming. It was published
by the Danish magazine Alt om DATA.

granular1d.c
------------
Simulation of granular media (Knudsen gas) in one dimension.
Test: `granular1d -T 1.0 -t 1.0 -N 100 -n 1000 -d 0.8 -R 0.01 -e 1.0`

License
-------
The package is distributed under GNU General Public License v2. See the file LICENSE for more
information.
