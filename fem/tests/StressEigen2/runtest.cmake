include(test_macros)
# No ElmerGrid step here on purpose. CMakeLists.txt copies the ready-made
# "angle" mesh, including the partitioning.4 and partitioning.8 directories
# that the NPROCS 4 and 8 runs need, and regenerating the mesh would leave
# those partitions describing the previous one.
#
# Running "ElmerGrid 1 2 angle.grd" here was also actively harmful, because
# angle.grd is not among the files CMakeLists.txt copies. Finding no input,
# ElmerGrid helpfully writes a built-in *3D* example in its place ("Because
# file angle.grd didn't exist, it was created for you"), which the next run
# of the test then meshes, replacing the 2D quad mesh this case is written
# for with a 3D hexahedral one. The first run passed and every run after it
# in the same build tree failed. CI never saw it because it always builds
# from scratch.
RUN_ELMER_TEST()
