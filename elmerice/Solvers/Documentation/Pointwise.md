# Pointwise
Pointwise.f90 was written as a way of interpolating scattered data on to the Elmer mesh. It may be superseded by [Grid2DInterpolator](./Grid2DInterpolator.md) (though it uses a different algorithm) now. It is not properly documented yet, but see comments in the pointwise.f90 file for information.
Note, if pointwise is used several times within one sif-file, the corresponding .so file needs to be duplicated under a different name (see also [Flowdepth](./FlowDepth.md) Solver), i.e. pointwise1.so and pointwise2.so.
