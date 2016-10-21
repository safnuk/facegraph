3D-Scan to Ribbon Graph Conversion
==================================

Given a triangulated surface with at least one boundary, and of negative Euler characteristic (i.e. if the surface has genus 0, then it must have at least three holes), this program
first performs combinatorial Ricci flow to find the conformally equivalent hyperbolic metric on the surface. Then, it uses the cut-locus construction to produce a metric ribbon graph embedded on the surface. Note that this ribbon graph allows one to recover the conformal structure (it is a *complete* conformal invariant of the surface). The primary application is for facial recognition - given a 3d scan of a face, the associated ribbon graph is a provides a (somewhat crude) signature that can be used for facial recognition. More importantly, it provides a method to canonically project the scan to a 2D-surface, which allows for more sophisticated algorithms to be run on the scan.

Installation and Compiling
--------------------------

The program depends on the [L-BFGS library](http://www.chokkan.org/software/liblbfgs/), which should be first installed on the system. It also uses [SparseLib++](http://math.nist.gov/sparselib++/) and [IML++](http://math.nist.gov/iml++/) libraries, but copies of these are included with the source code. To install and compile, run:
```
git clone https://github.com/safnuk/facegraph.git
cd facegraph/src
make main
```


Usage
-----

The program should run on any triangulated surface that has at least one boundary and negative euler characteristic. Run
```
main file.mesh
```
on such a file to produce the ribbon graph. Sample triangulated files are included in the data/ directory to get you started.

Demo
----

A [video](demo/cut-locus-demo.mp4) showing how the cut-locus construction works to produce a ribbon graph is in the demo/ directory.
