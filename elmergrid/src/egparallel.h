/* femelmer.h -> egparallel.h */
int FuseSolutionElmerPartitioned(char *prefix,char *outfile,int decimals,int parts,
				 int minstep, int maxstep, int dstep, int info);
int PartitionSimpleElements(struct FemType *data,struct ElmergridType *eg,struct BoundaryType *bound,
			    int dimpart[],int dimper[],int partorder, Real corder[],
			    Real parttol, int info);
int PartitionSimpleElementsNonRecursive(struct FemType *data,
					int dimpart[],int dimper[],int info);
#if USE_METIS
int PartitionConnectedElementsMetis(struct FemType *data,struct BoundaryType *bound,
				    struct ElmergridType *eg, int nparts,int metisopt,int info);
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
			 int partorder, Real corder[],Real parttol,int info);
#if USE_METIS
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
			      char *prefix,int decimals,int *parthalo,int indirect,
			      int parthypre,int subparts,int nooverwrite, int info);
