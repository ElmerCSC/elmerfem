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

typedef struct {
  int *Lcols, *Lrows, init;
  double *Lvals;
} ElmerRocalution;
ElmerRocalution S;


extern "C" void ROCSolve(int *n, int *rows, int *cols, double *vals, double *b, double *x_out, int *nonlin_update,int *fcomm )
{
using namespace rocalution;

    int i, *Lrows, *Lcols;
    double *Lvals;

    // Initialize rocALUTION
    init_rocalution();

    // rocALUTION objects
    LocalMatrix<double> mat;
    LocalVector<double> x, rhs;

    // Start time measurement
    double tick, tack;
    tick = rocalution_time();

    Lrows = (int *)malloc((*n+1)*sizeof(int));
    Lcols = (int *)malloc((rows[*n])*sizeof(int));
    Lvals= (double *)malloc((rows[*n])*sizeof(double));

    for(i=0; i<=*n; i++ ) Lrows[i] = rows[i];
    for(i=0; i<rows[*n]; i++ ) Lcols[i] = cols[i];
    for(i=0; i<rows[*n]; i++ ) Lvals[i] = vals[i];
    mat.SetDataPtrCSR(&Lrows, &Lcols, &Lvals, "A", rows[*n], *n, *n);

    // Move objects to accelerator
    mat.MoveToAccelerator();
    x.MoveToAccelerator();
    rhs.MoveToAccelerator();

    // Allocate vectors
    x.Allocate("x", mat.GetN());
    rhs.Allocate("rhs", mat.GetM());

    rhs.CopyFromData(b);
    x.CopyFromData(x_out);

    // Linear Solver
    IterativeLinearSolver<LocalMatrix<double>, LocalVector<double>, double> *ls;

    CG<LocalMatrix<double>, LocalVector<double>, double> ls_cg;
    BiCGStab<LocalMatrix<double>, LocalVector<double>, double> ls_bcg;
    BiCGStabl<LocalMatrix<double>, LocalVector<double>, double> ls_bcgl;
    GMRES<LocalMatrix<double>, LocalVector<double>, double> ls_gmres;
    FGMRES<LocalMatrix<double>, LocalVector<double>, double> ls_fgmres;
    IDR<LocalMatrix<double>, LocalVector<double>, double> ls_idr;

    int iterator = 2;
    switch(iterator) {
      case(0): ls = &ls_cg;  break;
      case(1): ls = &ls_bcg; break;
      case(2): ls = &ls_bcgl; ls_bcgl.SetOrder(4); break;
      case(3): ls = &ls_gmres; ls_gmres.SetBasisSize(50); break;
      case(4): ls = &ls_fgmres; ls_fgmres.SetBasisSize(50); break;
      case(5): ls = &ls_idr; ls_idr.SetShadowSpace(4); break;
      default: ls = &ls_bcg; break;
    }


    {
      double TOL=1e-6, bnrm;
      int MaxIter = 1000;
      bnrm = 0;
      for(i=0; i<*n; i++ ) bnrm += b[i]*b[i];
      bnrm = sqrt(bnrm);
      if (bnrm<1e-16) bnrm = 1;
      ls->Init(TOL*bnrm,1e-20,1e20,MaxIter);
    }


    // Preconditioner
    int prec = 0;

    Jacobi<LocalMatrix<double>, LocalVector<double>, double> prec_j;
    ILU<LocalMatrix<double>, LocalVector<double>, double> prec_i;

    switch(prec) {
      case(0): ls->SetPreconditioner(prec_j); break;

      case(1):
        {
          int level=0;
          prec_i.Set(level);
        }
        ls->SetPreconditioner(prec_i);
      break;
    }

    // Set solver operator
    ls->SetOperator(mat);

    // Build solver
    ls->Build();

    // Verbosity output
    ls->Verbose(2);

    // Print matrix info
    mat.Info();

    // Solve A x = rhs
    ls->Solve(rhs, &x);

    x.CopyToData(x_out);

    // Stop time measurement
    tack = rocalution_time();
    std::cout << "Solver execution:" << (tack - tick) / 1e6 << " sec" << std::endl;

    // Clear solver
    ls->Clear();

    // Stop rocALUTION platform
    stop_rocalution();

    free(Lrows); free(Lcols); free(Lvals);
}
