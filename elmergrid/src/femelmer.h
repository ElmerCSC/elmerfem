/* femelmer.h */
/* Routines for input, output and manipulation of Funcs and ElmerPost
   formats (programs created by Juha Ruokolainen at CSC). */ 
#define PARTMETIS 1

int LoadSolutionElmer(struct FemType *data,int results,char *prefix,int info);
int LoadElmerInput(struct FemType *data,struct BoundaryType *bound,
		   char *prefix,int nonames, int info);
int FuseSolutionElmerPartitioned(char *prefix,char *outfile,int decimals,int parts,
				 int minstep, int maxstep, int dstep, int info);
int SaveSolutionElmer(struct FemType *data,struct BoundaryType *bound,
		      int nobound,char *prefix,int decimals,int info);
int SaveSolutionElmerTriangles(struct FemType *data,char *prefix,int info);
int SaveElmerInput(struct FemType *data,struct BoundaryType *bound,
		   char *prefix,int decimals,int nooverwrite, int info);
int SaveSizeInfo(struct FemType *data,struct BoundaryType *bound,
		 char *prefix,int info);
int SaveElmerInputFemBem(struct FemType *data,struct BoundaryType *bound,
			 char *prefix,int decimals, int info);
int PartitionSimpleElements(struct FemType *data,struct ElmergridType *eg,struct BoundaryType *bound,
			    int dimpart[],int dimper[],int partorder, Real corder[],int info);
int PartitionSimpleElementsNonRecursive(struct FemType *data,
					int dimpart[],int dimper[],int info);
#if PARTMETIS
int PartitionConnectedElementsMetis(struct FemType *data,struct BoundaryType *bound,
				    int nparts,int metisopt,int info);
#endif
int ExtendBoundaryPartitioning(struct FemType *data,struct BoundaryType *bound,
			       int elemlayers,int info);
int PartitionSimpleElementsRotational(struct FemType *data,int dimpart[],int dimper[],
				      int info);
int PartitionConnectedElementsStraight(struct FemType *data,struct BoundaryType *bound,
				       struct ElmergridType *eg, int info);
int PartitionConnectedElements1D(struct FemType *data,struct BoundaryType *bound,
				 struct ElmergridType *eg, int info);
int PartitionSimpleNodes(struct FemType *data,int dimpart[],int dimper[],
			 int partorder, Real corder[],int info);
int LinearNodes(int elemtype);
#if PARTMETIS
int PartitionMetisMesh(struct FemType *data,struct ElmergridType *eg,
		       int partitions,int dual,int info);
int PartitionMetisGraph(struct FemType *data,struct BoundaryType *bound,
			struct ElmergridType *eg,int partitions,int metisopt,
			int dual,int info);
int ReorderElementsMetis(struct FemType *data,int info);
#endif
int OptimizePartitioningAtBoundary(struct FemType *data,struct BoundaryType *bound,int info);
int OptimizePartitioning(struct FemType *data,struct BoundaryType *bound,int noopt,
			 int partbw,int info);
int SaveElmerInputPartitioned(struct FemType *data,struct BoundaryType *bound,
			      char *prefix,int decimals,int halomode,int indirect,
			      int parthypre,int subparts,int nooverwrite, int info);
