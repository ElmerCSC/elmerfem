/* femelmer.h + femmesh.h -> egnative.h */
void InitGrid(struct GridType *grid);
void CreateExampleGrid(int dim,struct GridType **grids,int *nogrids,int info);
void SetElementDivision(struct GridType *grid,Real relh,int info);
void SetCellData(struct GridType *grid,struct CellType *cell,int info);
void CreateCells(struct GridType *grid,struct CellType **cell,int info);
void DestroyCells(struct CellType **cell); 
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
int LoadElmergrid(struct GridType **grid,int *nogrids,char *prefix,Real relh,int info);

void InitParameters(struct ElmergridType *eg);
int InlineParameters(struct ElmergridType *eg,int argc,char *argv[]);
int LoadCommands(char *prefix,struct ElmergridType *eg,
		 struct GridType *grid, int mode,int info);

int LoadElmerInput(struct FemType *data,struct BoundaryType *bound,
		   char *prefix,int nonames, int info);
int SaveElmerInput(struct FemType *data,struct BoundaryType *bound,
		   char *prefix,int decimals,int nooverwrite, int info);
int CreateElmerGridMesh(struct GridType *grid,
			struct FemType *data,struct BoundaryType *boundaries,
			Real relh,int info);
