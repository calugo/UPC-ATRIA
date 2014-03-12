
# README
This repository contains source codes to produce the results of the human atrial
cell model
reported in the article:
*Are SR Ca content fluctuations or SR refractoriness the key to atrial cardiac
alternans?: insights from a human atrial model (DOI:
10.1152/ajpheart.00515.2013)*

It contains four files:

1. This file.
2. UPCATRIALMODEL.tar.gz-C version of the program.
3. A version of the gnu scientific library.
4. UPC-ATRIALCELL.f-Fortran version of the program.

## The gnu scientific library:

An updated version of the  gsl library alongside installation instructions can
be obtained
at: http://www.gnu.org/software/gsl

Here it is only distributed for completeness, because the code heavily relies on
the library, so it should be installed in your system.

## The C version bundle: UPCATRIALMODEL.tar.gz

Once you have extracted the files in the compressed file you will find:

1. The source code itself, which is contained in the file UPCatrialCell.c
1. A folder called SAMPLEOUTPUT, which contains results obtained by running the
code as it is packed in the bundle.
1. A shell script in a file called 'atriacompile'

The parameters chosen to produce the sample results are:

* Stimulus period Ts=400 ms.
* Recovery times of 200ms (healthy cell) and 700ms for the alternating cell.

All the output files consist of three columns which correspond to time, healthy
cell and
alternating cell signals respectively.

## Compilation and Excecution (Unix/Mac OS).

The script 'atriacompile' executes the line:

`$cc -I/usr/local/include -o AtrialCell UPCatrialCell.c -L/usr/local/lib -lgsl
-lgslcblas -lm`

If the installation paths of the library are different they can be found by
executing:

`$gsl-config --libs --cflags`

If the code links and compiles correctly it generates the file:  `AtrialCell`,
which should
be run in a terminal as:

`$ ./AtrialCell`

## Fortran Code- UPC-ATRIALCELL.f

A fortran version of the code is also available which only needs to be compiled
and executed.
