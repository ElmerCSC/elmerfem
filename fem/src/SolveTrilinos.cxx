/* interface to Trilinos

This interface can in principal be used both in parallel and serially,
but currently it is only called in SParIterSolver.
It is compiled into Elmer by adding -DHAVE_TRILINOS
and linking with the Trilinos libraries. We currently
allow creating an iterative solver (from the Belos library)
and a preconditioner (ifpack or ML). To use a direct solver,
set "Ifpack Method" to "Amesos" and "Iterative Solver" to "None".

Older versions (below 10.10) of Trilinos demand a different call 
of the XML input file. This can be force by -DOLD_TRILINOS in the
CXXFLAGS 

To activate these linear solvers, set
'Linear System Use Trilinos' = Logical True
'Trilinos Input File' = String <xml filename>

see elmerfem/fem/examples/trilinos for an example.

*/

#include "../config.h"


#ifdef HAVE_TRILINOS

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef CHECK_ZERO
#undef CHECK_ZERO
#endif

// enable this to store matrices, generate debugging output etc
//#define DEBUG_TRILINOS_INTERFACE
//#define DUMP_IN_TRILINOS_INTERFACE

#define CHECK_ZERO(funcall) {ierr = funcall;\
if (ierr) {std::cout<<"Trilinos Error "<<ierr<<" returned from call "<<#funcall<<std::endl; return;}}

#define FCHECK_ZERO(funcall) {ierr = funcall;\
if (ierr) {std::cout<<"Trilinos Error "<<ierr<<" returned from call "<<#funcall<<std::endl; return Teuchos::null;}}

#define ERROR(msg,file,line)  {std::cerr << "Error in file "<<file<<", line "<<line<<":"<<std::endl; \
  std::cerr << msg << std::endl; \
  return;}

#define FERROR(msg,file,line)  {std::cerr << "Error in file "<<file<<", line "<<line<<":"<<std::endl; \
  std::cerr << msg << std::endl; \
  return Teuchos::null;}

#define WARNING(msg,file,line)  {std::cerr << "Warning in file "<<file<<", line "<<line<<":"<<std::endl; \
  std::cerr << msg << std::endl;}

#define PRINTLEVEL 0

#define OUT(s) if (am_printer) std::cout << s << std::endl;

#ifdef HAVE_MPI
#define HAD_MPI
#undef HAVE_MPI
#endif

#ifdef HAVE_HYPRE
#define HAD_HYPRE
#undef HAVE_HYPRE
#endif

#include "Epetra_SerialComm.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Export.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Operator.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Utils.hpp"

#include "BelosEpetraAdapter.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_StandardCatchMacros.hpp"

//#ifdef DEBUG_TRILINOS_INTERFACE
// only for debugging
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_BlockMapOut.h"
//#endif

#ifdef HAD_MPI
#define HAVE_MPI
#endif
#ifdef HAD_HYPRE
#define HAVE_HYPRE
#endif


typedef double ST;
typedef Teuchos::ScalarTraits<ST> SCT;
typedef SCT::magnitudeType MT;
typedef Epetra_MultiVector MV;
typedef Epetra_Operator OP;
typedef Belos::MultiVecTraits<ST,MV> MVT;
typedef Belos::OperatorTraits<ST,MV,OP> OPT;

#ifdef DEBUG_TRILINOS_INTERFACE
class MemTest 
  {
  public:
  
  MemTest() {std::cerr <<  "Elmer Trilinos container constructed"<<std::endl;}
  ~MemTest() {std::cerr << "Elmer Trilinos container destroyed"<<std::endl;}
  };
#endif

typedef struct ElmerTrilinosContainer {

Teuchos::RCP<Epetra_Comm> comm_;
Teuchos::RCP<Epetra_Map> assemblyMap_; // map with 'overlap' of nodes
Teuchos::RCP<Epetra_Map> solveMap_; // map with each node on one proc
Teuchos::RCP<Epetra_Export> exporter_;
Teuchos::RCP<Epetra_CrsMatrix> matrix_;
Teuchos::RCP<Epetra_Vector> rhs_;
Teuchos::RCP<Epetra_Vector> sol_;
Teuchos::RCP<Epetra_MultiVector> coords_; // node coordinates can be used by ML

double scaleFactor_; // scale entire system by a scalar constant, scaleFactor*Ax = scaleFactor*b
                     // enabled by setting "Scale Factor" in xml file

Teuchos::RCP<Teuchos::ParameterList> params_;

Teuchos::RCP<Epetra_Operator> prec_;
Teuchos::RCP<Belos::SolverManager<ST,MV,OP> > solver_;

Teuchos::RCP<struct ElmerTrilinosContainer> previous_;
Teuchos::RCP<struct ElmerTrilinosContainer> next_;

#ifdef DEBUG_TRILINOS_INTERFACE
Teuchos::RCP<MemTest> memtest_;
#endif
} ElmerTrilinosContainer;

// to avoid warnings from RCP's in debug mode when 
// temporarily returning the control to Fortran,   
// we store an extern pointer to a (doubly) linked 
// list of all 'Container' objects here
static Teuchos::RCP<ElmerTrilinosContainer> containerListHead=Teuchos::null;



// some auxiliary functions implemented below:

// creates the map without overlap
Teuchos::RCP<Epetra_Map> createSolveMap
        (Teuchos::RCP<Epetra_Comm> comm,
        int n, int* GID, int* owner);

// creates the matrix with overlap (e.g. shared nodes)
Teuchos::RCP<Epetra_CrsMatrix> createMatrix
        (Teuchos::RCP<Epetra_Map> assemblyMap,
        int *rows, int *cols, double* vals);

Teuchos::RCP<Epetra_Operator> createPreconditioner(
        Teuchos::RCP<Epetra_CrsMatrix> A,
        Teuchos::ParameterList& params,
        Teuchos::RCP<Epetra_MultiVector> coords);

Teuchos::RCP<Epetra_Operator> createIfpackPreconditioner(
        Teuchos::RCP<Epetra_CrsMatrix> A, Teuchos::ParameterList& params);

Teuchos::RCP<Epetra_Operator> createMLPreconditioner(
        Teuchos::RCP<Epetra_CrsMatrix> A, Teuchos::ParameterList& params,
        Teuchos::RCP<Epetra_MultiVector> coords);

Teuchos::RCP<Belos::SolverManager<ST,MV,OP> > createSolver
        (Teuchos::RCP<OP> A, Teuchos::RCP<OP> P,
        Teuchos::RCP<MV> x, Teuchos::RCP<MV> b,
        Teuchos::ParameterList& params);

static bool am_printer = true; // for output


// need to declare this extern so that the function names are
// consistent with Fortran
extern "C" {


// the calling order of SolveTrilinos1..4 is the same as for SolveHYPRE1..4

// construct matrix, solver and preconditioner
// nrows - number of local rows
// ncols - number of column indices on local partition
void SolveTrilinos1
 (
  int *n, int *nnz,
  int *rows, int *cols, double *vals,
  int *globaldof, int *owner, char* xmlfile, int* verbosityPtr, int** ContainerPtr,
  int *num_nodes, 
  double* xcoords, double* ycoords, double* zcoords,
  int *returnCode)
{

   bool serial=false;
   
   int verbose=*verbosityPtr;
   
   int& ierr=*returnCode;
   ierr=0;
   bool success=true;

   OUT("starting Trilinos setup");
   
   // Elmer has the unpleasant habit of finalizing MPI
   // if only one process is involved, this has to be 
   // changed in the source code if we want to use 
   // Trilinos for sequential runs as well. Check to make sure:
   int mpi_state;
   MPI_Initialized(&mpi_state);
   if (mpi_state==0)
     {
     OUT("MPI_Init has not been called, using serial comm");
     serial=true;
     }
   else
     {
     MPI_Finalized(&mpi_state);
     if (mpi_state!=0)
       {
       OUT("MPI_Finalize has already been called, using serial comm");
       serial=true;
       }
     }

   Teuchos::RCP<Epetra_Comm> comm;
   
   // create communicator
   if (serial)
     {
     comm = Teuchos::rcp(new Epetra_SerialComm());
     }
   else
     {
     comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
     }
   
   am_printer = (comm->MyPID()==0)&&(verbose>PRINTLEVEL);
   
   // sanity check, we expect standard CRS format
   if (*nnz != rows[*n])
     {
     ERROR("number of nonzeros incorrect",__FILE__,__LINE__);
     }

   Teuchos::RCP<Epetra_Map> assemblyMap;
   Teuchos::RCP<Epetra_CrsMatrix> A;

   // first create a map which has nodes shared between partitions
   // (e.g. an 'Elmer'-map)
   try {
   assemblyMap = Teuchos::rcp
        (new Epetra_Map(-1, *n,globaldof,1,*comm));
   } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success)

   if (!success)
     {
     WARNING("Failed to construct Elmer map",__FILE__,__LINE__);
     ierr = -1;
     return;
     }


#ifdef DUMP_IN_TRILINOS_INTERFACE 
CHECK_ZERO(EpetraExt::BlockMapToMatrixMarketFile("ElmerMap.txt",*assemblyMap));
#endif

   // based on this map, create a sparse matrix. 
   Teuchos::RCP<Epetra_CrsMatrix> A_elmer;
   try {
   A_elmer = 
        createMatrix(assemblyMap, rows, cols, vals);
   } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success)

   if (!success || A_elmer==Teuchos::null)
     {
     WARNING("Failed to construct Elmer matrix",__FILE__,__LINE__);
     ierr = -1;
     return;
     }

#ifdef DUMP_TRILINOS_INTERFACE
   CHECK_ZERO(EpetraExt::RowMatrixToMatrixMarketFile("ElmerMatrix.txt",*A_elmer));
#endif
     
   // now construct a map with each node owned by one partition (a solver-map)
   Teuchos::RCP<Epetra_Map> solveMap;
   try {   
    solveMap = createSolveMap(comm,*n,globaldof,owner);
   } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success)

   if (!success || solveMap==Teuchos::null)
     {
     WARNING("Failed to construct map",__FILE__,__LINE__);
     ierr = -1;
     return;
     }

#ifdef DUMP_IN_TRILINOS_INTERFACE 
CHECK_ZERO(EpetraExt::BlockMapToMatrixMarketFile("TrilinosMap.txt",*solveMap));
#endif

   Teuchos::RCP<Epetra_Export> exporter;

   // construct an exporter to transfer data between the two object types
   try {
   exporter = Teuchos::rcp(new Epetra_Export
        (*assemblyMap, *solveMap));
   } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success)

   if (!success || exporter==Teuchos::null)
     {
     WARNING("Failed to construct exporter",__FILE__,__LINE__);
     ierr = -1;
     return;
     }


try {   
   // build the non-overlapping matrix
   A = Teuchos::rcp(new Epetra_CrsMatrix
        (Copy, *solveMap, A_elmer->MaxNumEntries(),false));
   // create the matrix from the overlapping input matrix
   CHECK_ZERO(A->Export(*A_elmer, *exporter, Add));
   CHECK_ZERO(A->FillComplete());

   } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success)
   
   if (!success)
     {
     WARNING("Failed to construct matrix",__FILE__,__LINE__);
     ierr = -1;
     return;
     }
   
   // Create the rhs and solution
   Teuchos::RCP<Epetra_Vector> sol=Teuchos::rcp(new Epetra_Vector(A->RowMap()));
   Teuchos::RCP<Epetra_Vector> rhs=Teuchos::rcp(new Epetra_Vector(A->RowMap()));

   OUT("matrix constructed")


//////////////////////////////////////////////////////////////
// read parameters from input file 
//////////////////////////////////////////////////////////////

   // read solver parameters from a file
   Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList());
   std::string filename(xmlfile);
   if (filename=="none")
     {
     WARNING("no parameter file specified, using default settings in Trilinos",__FILE__,__LINE__);
     }
   else
     {
     OUT("reading parameters from '"+filename+"'");
     try {
#ifdef OLD_TRILINOS
       Teuchos::updateParametersFromXmlFile(filename,params.get());
#else
// for Trilinos 10.10 and later     
       Teuchos::updateParametersFromXmlFile(filename,params.ptr());
#endif       
       } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);
       if (!success)
         {
         WARNING("failed to read your Trilinos input file, using default values",
                __FILE__,__LINE__);
         ierr = 1;
         }
     }
     
   bool print_matrix = params->get("Dump Matrix",false);

   if (print_matrix)
   {
#if 1
//#ifdef DEBUG_TRILINOS_INTERFACE
   std::string filename = params->get("Filename Base","Trilinos")+"Matrix.mtx";
   CHECK_ZERO(EpetraExt::RowMatrixToMatrixMarketFile(filename.c_str(),*A));
#else
   WARNING("you have specified 'Dump Matrix', but DEBUG_TRILINOS_INTERFACE is not defined",
        __FILE__,__LINE__);
#endif
   }

   double scaleFactor = params->get("Scale Factor",1.0);
   if (scaleFactor!=1.0) CHECK_ZERO(A->Scale(scaleFactor));

   Teuchos::RCP<Epetra_MultiVector> coords0=Teuchos::rcp(new Epetra_MultiVector(*assemblyMap,3));
   Teuchos::RCP<Epetra_MultiVector> coords=Teuchos::rcp(new Epetra_MultiVector(*solveMap,3));
   
   int k = *num_nodes;
   int dof = (int)(assemblyMap->NumMyElements()/k);
   if (dof*k != assemblyMap->NumMyElements())
     {
     ERROR("size mismatch of coord arrays",__FILE__,__LINE__);
     }

   for (int i=0;i<k; i++)
     {
     for (int j=0;j<dof;j++)
       {
       (*coords0)[0][dof*i+j] = xcoords[i];
       (*coords0)[1][dof*i+j] = ycoords[i];
       (*coords0)[2][dof*i+j] = zcoords[i];
       }
     }

   CHECK_ZERO(coords->Export(*coords0, *exporter, Zero));
   
   Epetra_Map auxMap(-1,k,0,*comm);
   Teuchos::RCP<Epetra_MultiVector> auxVec = Teuchos::rcp(new Epetra_MultiVector(auxMap,3));
   for (int i=0;i<k;i++)
     { 
     for (int j=0;j<3;j++)
       {
       (*auxVec)[j][i] = (*coords)[j][dof*i];
       }
     }
   
   coords = auxVec;
   
  ///////////////////////////////////////////////////////////////////
  // create/setup preconditioner                                   //
  ///////////////////////////////////////////////////////////////////
  Teuchos::RCP<Epetra_Operator> prec = createPreconditioner(A, *params, coords);

  //////////////////////////////////////////////////////////////
  // Krylov subspace method setup (Belos)                     //
  //////////////////////////////////////////////////////////////
  Teuchos::RCP<Belos::SolverManager<ST,MV,OP> > solver = createSolver(A,prec,
        sol, rhs, *params);

  // print parameter list so we can see default values
  // and if some of our input is unused.
  if (am_printer && verbose>=5)
    {
    std::cout << "Trilinos parameter list: "<<std::endl;
    std::cout << *params<<std::endl;
    }


   //////////////////////////////////////////////////////////////
   // construct a container to return to Fortran               //
   //////////////////////////////////////////////////////////////
   if (*ContainerPtr != 0)
     {
     WARNING("pointer passed into SolveTrilinos1 not NULL, possible memory leak.",
        __FILE__,__LINE__);
     }

   Teuchos::RCP<ElmerTrilinosContainer> Container = 
        Teuchos::rcp(new ElmerTrilinosContainer());

   *ContainerPtr=(int*)(Container.get());
   
   // store a pointer to the container in a static 
   // linked list to avoid confusing the garbage   
   // collection (Teuchos RCPs won't understand that
   // we have stored a raw pointer in Elmer somewhere)
   if (containerListHead!=Teuchos::null)
      {
      containerListHead->previous_=Container;
      }
    Container->next_=containerListHead;
    containerListHead = Container;
   
   // put pointers to all the important objects 
   // in the container.
   Container->comm_=comm;
   Container->params_=params;
   Container->assemblyMap_=assemblyMap;
   Container->solveMap_=solveMap;
   Container->exporter_=exporter;
   Container->matrix_=A;
   Container->coords_=coords;
   Container->scaleFactor_=scaleFactor;
   Container->rhs_=rhs;
   Container->sol_=sol;
  Container->prec_=prec;
  Container->solver_=solver;
  
#ifdef DEBUG_TRILINOS_INTERFACE
  Container->memtest_=Teuchos::rcp(new MemTest());
#endif
  
  return;

  }
////////////////////////////////////////////////////////////

void SolveTrilinos2
 (
  int *n, double *xvec, double *rhsvec, int *Rounds, double *TOL,
  int *verbosityPtr, int** ContainerPtr,
  int* returnCode
 )
  {
  int verbose = *verbosityPtr;
   int& ierr=*returnCode;   
   ierr=0;
   bool success=true;

ElmerTrilinosContainer* Container = (ElmerTrilinosContainer*)(*ContainerPtr);
if (Container==NULL) ERROR("invalid pointer passed to SolveTrilinos2",__FILE__,__LINE__);

   am_printer = (Container->comm_->MyPID()==0)&&(verbose>PRINTLEVEL);
   
   // get the data structures built in SolveTrilinos1():
   Teuchos::RCP<Epetra_Map> assemblyMap = Container->assemblyMap_;
   Teuchos::RCP<Epetra_Map> solveMap = Container->solveMap_;
   Teuchos::RCP<Epetra_Export> exporter = Container->exporter_;
   Teuchos::RCP<Epetra_CrsMatrix> A = Container->matrix_;
   Teuchos::RCP<Epetra_Vector> x = Container->sol_;
   Teuchos::RCP<Epetra_Vector> b = Container->rhs_;
   Teuchos::RCP<Epetra_Operator> prec = Container->prec_;
   Teuchos::RCP<Belos::SolverManager<ST,MV,OP> > solver = Container->solver_;
   Teuchos::RCP<Teuchos::ParameterList> params = Container->params_;
   double scaleFactor = Container->scaleFactor_;

   // import the vectors
   Epetra_Vector bview(View, *assemblyMap, rhsvec);
   CHECK_ZERO(b->Export(bview, *exporter, Add));

   // import the vectors
   Epetra_Vector xview(View, *assemblyMap, xvec);
   CHECK_ZERO(x->Export(xview, *exporter, Zero));

   if (scaleFactor!=1.0)
     { 
     CHECK_ZERO(b->Scale(scaleFactor));
     }

   // override the settings for tconv tol and num iter using Elmer inut data:
   if (*TOL>=0.0) params->sublist("Belos").set("Convergence Tolerance",*TOL);
   if (*Rounds>0) params->sublist("Belos").set("Maximum Iterations",*Rounds);
   
  // check initial residual - do not start solver if already converged.
  
  int numrhs=1;
  std::vector<double> actual_resids( numrhs );
  std::vector<double> rhs_norm( numrhs );
  Epetra_MultiVector resid(A->OperatorRangeMap(), numrhs);

#if 0
  OPT::Apply( *A, *x, resid );
  MVT::MvAddMv( -1.0, resid, 1.0, *b, resid );
  MVT::MvNorm( resid, actual_resids );
  MVT::MvNorm( *b, rhs_norm );  
  bool converged = true;
  
  for ( int i=0; i<numrhs; i++)
    {
    double actRes = actual_resids[i]/rhs_norm[i];
    if (actRes >= *TOL) converged = false;
    } 
  
  if (converged)
    {
    OUT("initial solution passed to Trilinos is already good enough - returning to Elmer.");
    return;
    }
#else
// always start with a 0 initiaial guess because Belos
// does a relative residual check only
//x->PutScalar(0.0);
#endif   

/////////////////////////////////////////////////////////////
// Perform solve                                           //
/////////////////////////////////////////////////////////////

  if (solver!=Teuchos::null)
    {
    OUT("start iterative solver");
    solver->setParameters(Teuchos::rcp(&(params->sublist("Belos")),false));
    solver->reset(Belos::Problem);
    Belos::ReturnType ret;
    try {
    ret = solver->solve();
    } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);

if (!success) 
  {
  WARNING("Trilinos solve failed!",__FILE__,__LINE__);
  ierr=-1; // caught an exception -> error
  return;
  }

    // check for loss of accuracy
    bool loa = solver->isLOADetected();

    if (loa)
      {
      WARNING("loss of accuracy in Belos solve!",__FILE__,__LINE__);
      ierr=1;
      }

    if (ret!=Belos::Converged)
      {
      WARNING("Belos did not converge!",__FILE__,__LINE__);
      ierr=2;
      }
    }
  else if (prec!=Teuchos::null)
    {
    OUT("apply operator inverse");
    CHECK_ZERO(prec->ApplyInverse(*b,*x));
    }
  else
    {
   WARNING("no solver or preconditioner available",__FILE__,__LINE__);
    ierr=3;
    *x=*b;
    }

   bool print_vectors = params->get("Dump Vectors",false);
   if (print_vectors)
   {
#ifdef DEBUG_TRILINOS_INTERFACE   
   std::string filebase = params->get("Filename Base","Trilinos");
   EpetraExt::MultiVectorToMatrixMarketFile((filebase+"Rhs.txt").c_str(),*b);
   EpetraExt::MultiVectorToMatrixMarketFile((filebase+"Sol.txt").c_str(),*x);
#else
   WARNING("you have specified 'Dump Vectors', but DEBUG_TRILINOS_INTERFACE is not defined",
        __FILE__,__LINE__);
#endif
   }

   

#ifdef DEBUG_TRILINOS_INTERFACE   
    
  //
  // Compute actual residual
  //
  bool badRes = false;

  OPT::Apply( *A, *x, resid );
  MVT::MvAddMv( -1.0, resid, 1.0, *b, resid );
  MVT::MvNorm( resid, actual_resids );
  MVT::MvNorm( *b, rhs_norm );  

  if (verbose>=3)
    {
    std::cout << "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
    }
  for ( int i=0; i<numrhs; i++) 
    {
    double actRes = actual_resids[i]/rhs_norm[i];
    if (verbose>=3)
      {
      std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
      }
    if (actRes > *TOL) badRes = true;
    } 

  if (badRes) 
    {
    WARNING("bad actual residual found!",__FILE__,__LINE__);
    std::cerr << "(required accuracy was "<<*TOL<<")"<<std::endl;
    ierr=4;
    }

#endif

   // import the vectors
   CHECK_ZERO(xview.Import(*x, *exporter, Zero));

// Trilinos cleans up itself (because of RCP's)
return;
}

////////////////////////////////////////////////////////////
// destructor                                             //
////////////////////////////////////////////////////////////

void SolveTrilinos4(int** ContainerPtr)
  {
  ElmerTrilinosContainer* Container = 
        (ElmerTrilinosContainer*)(*ContainerPtr);

#ifdef DEBUG_TRILINOS_INTERFACE  
std::cerr << "PID "<<Container->comm_->MyPID()<<": destroy Trilinos object "<<std::endl;
#endif

  // remove this container from the list
  if (Container->next_!=Teuchos::null)
    {
    Container->next_->previous_=Container->previous_;
    }
  if (Container->previous_!=Teuchos::null)
    {
    Container->previous_->next_=Container->next_;
    }

  // nullify the pointer in fortran
  if (Container == containerListHead.get())
    {
    containerListHead=Container->next_;
    }

  *ContainerPtr=0;
  Container = NULL;
  
  
  
  // do some global barriers after destruction
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  }                                                                                                                                         

// this function deletes ALL Trilinos objects created by Elmer that     
// have not been destroyed properly, yet (by SolveTrilinos4).           
// It should only be called at the very end of an Elmer run.            
void TrilinosCleanup(void)
  {
#ifdef DEBUG_TRILINOS_INTERFACE  
  std::cout << "Destroying all remaining Trilinos objects..."<<std::endl;
#endif
  Teuchos::RCP<struct ElmerTrilinosContainer> c = containerListHead;
  while (c!=Teuchos::null)
    {
    c=c->next_;
    if (c!=Teuchos::null) c->previous_ = Teuchos::null;
    }
  containerListHead = Teuchos::null;
  }


}//extern "C"

// creates the map without overlap
Teuchos::RCP<Epetra_Map> createSolveMap
        (Teuchos::RCP<Epetra_Comm> comm,
        int n, int* GID, int* owner)
  {
  int ierr;
  // first figure out how many nodes are 'owned' by this process
  int nloc = 0;
  for (int i=0;i<n;i++) nloc+=owner[i];
  
  // construct an array with only owned global node id's
  int *MyGIDs = new int[nloc];
  int k=0;
  for (int i=0;i<n;i++)
    {
    if (owner[i])
      {
      MyGIDs[k++] = GID[i];
      }
    }
  
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp
        (new Epetra_Map(-1, nloc, MyGIDs, 1, *comm));
  
  delete [] MyGIDs;
  return map;
  }

// creates the matrix with overlap (e.g. shared nodes)
Teuchos::RCP<Epetra_CrsMatrix> createMatrix
        (Teuchos::RCP<Epetra_Map> assemblyMap,
        int *rows, int *cols, double* vals)
  {
  int ierr=0;
  bool success=true;
  Teuchos::RCP<Epetra_CrsMatrix> A;
  
  if (assemblyMap==Teuchos::null)
    {
    FERROR("map passed to createMatrix is null",__FILE__,__LINE__);
    }
  
  int nrows = assemblyMap->NumMyElements();
  
  int *row_size = new int[nrows];
  
  int max_row_size=0;
  
  for (int i=0;i<nrows;i++)
    {
    row_size[i]=rows[i+1]-rows[i];
    if (row_size[i]>max_row_size) max_row_size=row_size[i];
    }
  
  int *gcols=new int[max_row_size];

  A=Teuchos::rcp(new Epetra_CrsMatrix(Copy,*assemblyMap,row_size,true));
  
  // Now go through my local rows and set the matrix entries.
  try {
  for (int i = 0; i < nrows; i++)
    {
    for (int j=0;j<row_size[i];j++) gcols[j]=assemblyMap->GID(cols[rows[i]-1+j]-1);
    FCHECK_ZERO(A->InsertGlobalValues(assemblyMap->GID(i),row_size[i],&(vals[rows[i]-1]),gcols));  
    }
  } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success)
  delete [] row_size;
  delete [] gcols;
  
  // Assemble after setting the coefficients
  FCHECK_ZERO(A->FillComplete());
  return A;
  }


Teuchos::RCP<Epetra_Operator> createPreconditioner(
        Teuchos::RCP<Epetra_CrsMatrix> A, 
        Teuchos::ParameterList& params,
        Teuchos::RCP<Epetra_MultiVector> coords)
        {
        std::string type="None";
        type=params.get("Preconditioner",type);
        if (type=="Ifpack") return createIfpackPreconditioner(A,params);
        if (type=="ML") return createMLPreconditioner(A,params,coords);
        if (type=="None") return Teuchos::null;
        WARNING("invalid 'Preconditioner', returning null",__FILE__,__LINE__);
        return Teuchos::null;
        }


Teuchos::RCP<Epetra_Operator> createIfpackPreconditioner(
        Teuchos::RCP<Epetra_CrsMatrix> A,
        Teuchos::ParameterList& params)
  {
  int ierr;
  Teuchos::RCP<Epetra_Operator> prec = Teuchos::null;
  std::string type = "None";
  type=params.get("Preconditioner",type);
  
  if (type=="Ifpack")
    {
    Teuchos::RCP<Ifpack_Preconditioner> ifPrec;
    Teuchos::ParameterList& ifpackList = params.sublist("Ifpack");
    std::string ifpackType = "Amesos";
    ifpackType=params.get("Ifpack Preconditioner",ifpackType);
    OUT(ifpackType);
    int OverlapLevel = params.get("Ifpack Overlap Level",0); 
    
    OUT("construct ifpack preconditioner object");
    Ifpack factory;
    
    ifPrec = Teuchos::rcp(factory.Create(ifpackType,
               A.get(), OverlapLevel) );

    OUT("set parameters");
    ifPrec->SetParameters(ifpackList);
    OUT("initialize");
    FCHECK_ZERO(ifPrec->Initialize());
      
    OUT("compute");
    FCHECK_ZERO(ifPrec->Compute());                                                                                                 
    
    prec = ifPrec;
    }
  else
    {
    FERROR(" 'Preconditioner' must be 'Ifpack' for this function.",
        __FILE__,__LINE__);
    }
  return prec;
  }

Teuchos::RCP<Epetra_Operator> createMLPreconditioner(
        Teuchos::RCP<Epetra_CrsMatrix> A,
        Teuchos::ParameterList& params,
        Teuchos::RCP<Epetra_MultiVector> coords)
  {
  Teuchos::ParameterList& mlParams = params.sublist("ML");
 
  std::string defaults = "SA";
  if (mlParams.isParameter("default values"))
    {
    defaults = mlParams.get("default values",defaults);
    }
  else if (mlParams.isParameter("SetDefaults"))
    {
    defaults = mlParams.get("SetDefaults",defaults);
    }
  
  ML_Epetra::SetDefaults(defaults,mlParams,NULL,NULL,false);
  
  if (mlParams.get("aggregation: aux: enable",false)==true)
    {
    mlParams.set("x-coordinates",(*coords)[0]);
    mlParams.set("y-coordinates",(*coords)[1]);
    mlParams.set("z-coordinates",(*coords)[2]);
    }
  
  Teuchos::RCP<Epetra_Operator> prec = Teuchos::null;
  
  const Epetra_MpiComm *mpiComm = 
        dynamic_cast<const Epetra_MpiComm*>(&(A->Comm()));

  if (mpiComm==NULL) FERROR("ML needs an MpiComm object, run with ElmerSolver_mpi",__FILE__,__LINE__);

  OUT("create ML Preconditioner...");
  
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPreconditionerPtr = 
        Teuchos::rcp( new ML_Epetra::MultiLevelPreconditioner(*A, mlParams, true) );
  
  prec = MLPreconditionerPtr;                                       

  bool analyze_cycle = params.get("ML: Analyze Cycle",false);
  bool dump_matrices = params.get("ML: Dump Matrices",false);      
  bool test_smoothers = params.get("ML: Test Smoothers",false);

  if (test_smoothers || analyze_cycle)
    {
    // ML dumps all of its output to std::cout
    std::ostream& os = std::cout;
    int NumSmoo = mlParams.get("smoother: sweeps",1);
    int NumCycles = mlParams.get("cycle applications",1);
    int NumPre=0, NumPost=0;
    std::string PreOrPost = mlParams.get("smoother: pre or post","both");
    if ((PreOrPost=="pre")||(PreOrPost=="both"))
      {
      NumPre=NumCycles;
      }
    if ((PreOrPost=="post")||(PreOrPost=="both"))
      {
      NumPost=NumCycles;
      }
      
      
    if (analyze_cycle) 
      {
      if (prec->Comm().MyPID()==0)
        {
        os << "Print hierarchy and analyze effect of MG cycle on a random vector\n";
        }
      MLPreconditionerPtr->AnalyzeHierarchy(true,NumPre,NumPost,NumCycles);
      }
    if (test_smoothers)    
      {
      if (prec->Comm().MyPID()==0)
        {
        os << "Test various Smoothers"<<std::endl;
        }
      Teuchos::ParameterList testList(mlParams);
      testList.set("ML validate parameter list",false);
      testList.set("cycle applications",1); // used in a Kcylov method here
      testList.set("test: Jacobi",true);
      testList.set("test: Gauss-Seidel",true);
      testList.set("test: symmetric Gauss-Seidel",true);
      testList.set("test: block Gauss-Seidel",true);
      testList.set("test: Aztec",false); // causes an exception
      testList.set("test: IFPACK",true);
      testList.set("test: ParaSails",true);
      testList.set("test: ML self smoother",true);
      
      MLPreconditionerPtr->TestSmoothers(testList);
      }
    }
  if (dump_matrices) 
    {
    WARNING("'ML: Dump Matrices' not implemented",__FILE__,__LINE__);
    //DumpMLHierarchy(MLPreconditionerPtr);
    }
  return prec;
  }
                                                                                  
// create an iterative solver for the linear system Ax=b,
// preconditioned by P.
Teuchos::RCP<Belos::SolverManager<ST,MV,OP> > createSolver
        (Teuchos::RCP<OP> A, Teuchos::RCP<OP> P,
         Teuchos::RCP<MV> x, Teuchos::RCP<MV> b,
         Teuchos::ParameterList& params)
  {
  Teuchos::RCP<Belos::SolverManager<ST,MV,OP> > belosSolverPtr;


  // retrieve User's Belos list, add more things
  Teuchos::ParameterList& belosList=params.sublist("Belos");

  //! Belos preconditioner interface
  Teuchos::RCP<Belos::EpetraPrecOp> belosPrecPtr=Teuchos::null;

  //! Belos linear problem interface
  Teuchos::RCP<Belos::LinearProblem<ST,MV,OP> > belosProblemPtr;

  std::string linearSolver=params.get("Iterative Solver","GMRES");
  if (linearSolver=="None")
    {
    return Teuchos::null;
    }
  bool verbose = true;
#ifdef DEBUG_TRILINOS_INTERFACE
  bool debug = true;
#else
  bool debug = false;
#endif  

  int verbosity = Belos::Errors + Belos::Warnings;
  if (verbose)
    { //TODO: where to put which option? how do we get readable output?
    verbosity+=Belos::TimingDetails+Belos::IterationDetails;
    verbosity+=Belos::StatusTestDetails+Belos::OrthoDetails+Belos::FinalSummary;
    }
  if (debug) verbosity+=Belos::Debug;
  // User is allowed to override these settings
  if (belosList.isParameter("Verbosity")==false)
    {
    belosList.set("Verbosity",verbosity);
    }

  if (belosList.isParameter("Output Stream")==false)
    {
    belosList.set("Output Stream",Teuchos::rcp(&std::cout, false));
    }

  belosList.set("Output Style",(int)Belos::Brief);
  
  // by default, Belos checks only ||r_imp||_2/||r_0||_2, but if we're
  // solving a sequence of systems and therefore have a very good initial
  // guess, it is better to check ||r||_2/||b||_2 (actually a more fancy 
  // convergence test like the one used in Elmer would be better, but I  
  // haven't figured out how to tell Belos to do that)
  if (linearSolver=="GMRES" && belosList.isParameter("Implicit Residual Scaling")==false)
    {
    belosList.set("Implicit Residual Scaling","Norm of RHS");
    }

  // create Belos interface to preconditioner.
  // This is simply an Epetra_Operator with 'Apply' and 'ApplyInverse' switched.
  if (!Teuchos::is_null(P))
    {
    belosPrecPtr = Teuchos::rcp(new Belos::EpetraPrecOp(P));
    }
  // create Belos problem interface
  belosProblemPtr = Teuchos::rcp(new Belos::LinearProblem<ST,MV,OP>(A,x,b));

  if (belosPrecPtr!=Teuchos::null)
    {
    // set preconditioner
    OUT("set preconditioner...");
    belosProblemPtr->setLeftPrec(belosPrecPtr);
    }
  bool set = belosProblemPtr->setProblem();
  if (set == false) {
    FERROR("Belos::LinearProblem failed to set up correctly!",__FILE__,__LINE__);
    }

  // create the solver
  Teuchos::RCP<Teuchos::ParameterList> belosListPtr=rcp(&belosList,false);
  if (linearSolver=="CG")
    {
    belosSolverPtr = Teuchos::rcp(new Belos::BlockCGSolMgr<ST,MV,OP>(belosProblemPtr,belosListPtr));
    }
  else if (linearSolver=="GMRES")
    {
    Teuchos::RCP<Teuchos::ParameterList> belosListPtr=Teuchos::rcp(&belosList,false);
    belosSolverPtr = Teuchos::rcp(new Belos::BlockGmresSolMgr<ST,MV,OP>(belosProblemPtr,belosListPtr));
    }
  else if (linearSolver=="None")
    {
    belosSolverPtr = Teuchos::null;
    }
  else
    {
    FERROR("Currently 'CG', 'GMRES' and 'None' \n"
        " are supported as 'Iterative Solver'",__FILE__,__LINE__);
    }
  return belosSolverPtr;
  }

void TrilinosSetNodeCoords(int* num_nodes, double* x, double* y, double* z)
  {
  }

#endif
