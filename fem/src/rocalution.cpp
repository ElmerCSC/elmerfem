/* ************************************************************************
 * Copyright (C) 2018-2021 Advanced Micro Devices, Inc. All rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * ************************************************************************ */

#include <cstdlib>
#include <iostream>
#include <rocalution/rocalution.hpp>

extern "C" void ROCSolve(int *n, int *rows, int *cols, double *vals, double *b, double *x_out, int *nonlin_update,int *fcomm )
{
using namespace rocalution;

    int i;

    // Initialize rocALUTION
    init_rocalution();

    // Print rocALUTION info
//    info_rocalution();

    // rocALUTION objects
    LocalVector<double> x;
    LocalVector<double> rhs;
    LocalVector<double> e;
    LocalMatrix<double> mat;

    int *Lrows = new int[*n+1];
    int *Lcols = new int[rows[*n]];
    double *Lvals = new double[rows[*n]];

    for(i=0; i<=*n; i++ ) Lrows[i]=rows[i];
    for(i=0; i<rows[*n]; i++ ) Lcols[i]=cols[i];
    for(i=0; i<rows[*n]; i++ ) Lvals[i]=vals[i];

    // Read matrix from MTX file
//    mat.ReadFileMTX("xx");
    mat.SetDataPtrCSR(&Lrows, &Lcols, &Lvals, "A", rows[*n], *n, *n);


    // Move objects to accelerator
    mat.MoveToAccelerator();
    x.MoveToAccelerator();
    rhs.MoveToAccelerator();
    e.MoveToAccelerator();

    // Allocate vectors
    x.Allocate("x", mat.GetN());
    rhs.Allocate("rhs", mat.GetM());
    e.Allocate("e", mat.GetN());

    // Initialize rhs & x
//     for(i=0; i<*n; i++) { rhs[i] = b[i]; x[i] = x_out[i]; }

    rhs.CopyFromData(b);
    x.CopyFromData(x_out);

    // Linear Solver
    CG<LocalMatrix<double>, LocalVector<double>, double> ls;

    // Preconditioner
    Jacobi<LocalMatrix<double>, LocalVector<double>, double> p;


    // Initial zero guess
//    x.Zeros();

    // Set solver operator
    ls.SetOperator(mat);

    // Set solver preconditioner
    ls.SetPreconditioner(p);

    // Build solver
    ls.Build();

    // Verbosity output
    ls.Verbose(1000);

    // Print matrix info
    mat.Info();

    // Start time measurement
    double tick, tack;
    tick = rocalution_time();

    // Solve A x = rhs
    ls.Solve(rhs, &x);

//      LocalVector<double> x0;
//      x0.Allocate("x0", mat.GetN());
//      x0.CopyFrom(x);
      for(i=0; i<*n; i++) x_out[i] = x[i];

    // Stop time measurement
    tack = rocalution_time();
    std::cout << "Solver execution:" << (tack - tick) / 1e6 << " sec" << std::endl;

    // Clear solver
     ls.Clear();

   // Stop rocALUTION platform
     stop_rocalution();
}
