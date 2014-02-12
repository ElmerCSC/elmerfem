/* egnative.h */
/* Subroutines for creating and manipulating the native format of ElmerGrid. */

void InitGrid(struct GridType *grid);
void CreateExampleGrid(int dim,struct GridType **grids,int *nogrids,int info);
void SetElementDivision(struct GridType *grid,Real relh,int info);
void SetCellData(struct GridType *grid,struct CellType *cell,int info);
int SetCellKnots(struct GridType *grid, struct CellType *cell,int info);
int SetCellKnots1D(struct GridType *grid, struct CellType *cell,int info);
int GetKnotCoordinate(struct CellType *cell,int i,int j,Real *x,Real *y);
int GetKnotIndex(struct CellType *cell,int i,int j);
int GetElementIndices(struct CellType *cell,int i,int j,int *ind);
int GetElementIndex(struct CellType *cell,int i,int j);
int GetElementCoordinates(struct CellType *cell,int i,int j,
			  Real *globalcoord,int *ind);
int GetSideInfo(struct CellType *cell,int cellno,int side,int element,
		int *elemind);
void SetElementDivisionExtruded(struct GridType *grid,int info);
void SetElementDivisionCylinder(struct GridType *grid,int info);

int SaveElmergrid(struct GridType *grid,int nogrids,char *prefix,int info);
int LoadElmergrid(struct GridType **grid,int *nogrids,char *prefix,int info);
void InitParameters(struct ElmergridType *eg);
int LoadCommands(char *prefix,struct ElmergridType *eg,
		 struct GridType *grid, int mode,const char *IOmethods[],
		 int info);
int CreateElmerGridMesh(struct GridType *grid,
			struct FemType *data,struct BoundaryType *boundaries,
			Real relh,int info);
