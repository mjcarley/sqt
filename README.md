SQT is a simple library implementing (singular) quadrature methods, in
particular adaptive Gaussian quadrature, for three-dimensional
boundary integral methods. It supplies functions intended to be used
inside BEM codes, so some of the interfaces are quite basic but the
examples in ./tests/ show how to use the codes in practice.

# Prerequisites

You will need to have installed the blaswrap wrappers for BLAS

- https://github.com/mjcarley/blaswrap


# Installation

Installation instructions are given in the file INSTALL. If you have
downloaded the github distribution, you will need to set up the
automake system by running

./autogen.sh

in the root directory of the source code. You can then go through the
usual autotools configure process

./configure (options)

If you want to use code optimized with AVX extensions (this is
experimental but does give a performance improvement in the code where
it is in use) pass the compiler flag SQT_USE_AVX. For example:

CFLAGS="-O3 -g -DSQT_USE_AVX" ./configure ...

To compile the code:

make

and to install it:

make install

If you want to generate the documentation:

make doxygen-doc

installs the html documentation in doc/html
