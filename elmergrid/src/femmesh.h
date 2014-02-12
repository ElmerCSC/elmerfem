/* femmesh.h */
/* Subroutines for creating a mesh that operate on 
   variables that are of type GridType and CellType. 
   All subroutines operate on the sub-cell level - not on 
   individual knots or elements. */

void InitGrid(struct GridType *grid);
void ExampleGrid1D(struct GridType **grids,int *nogrids,int info);
void ExampleGrid2D(struct GridType **grids,int *nogrids,int info);
void ExampleGrid3D(struct GridType **grids,int *nogrids,int info);
void SetElementDivision(struct GridType *grid,Real relh,int info);
void SetCellData(struct GridType *grid,struct CellType *cell,int info);
void CreateCells(struct GridType *grid,struct CellType **cell,int info);
void DestroyCells(struct CellType **cell); 
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
int InspectElement(struct FemType *data,int idx);
