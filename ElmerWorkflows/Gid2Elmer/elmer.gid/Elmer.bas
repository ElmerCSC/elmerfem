==================================================================
                        Elmer Mesh Files
==================================================================

*set elems(all)
*loop nodes
*format "%8i -1 %16.8e %16.8e %16.8e"
mesh.nodes: *NodesNum *NodesCoord(1,real) *NodesCoord(2,real) *NodesCoord(3,real)
*end nodes

*set var Boundary202=0
*set var Boundary203=0
*set elems(Linear)
*set Cond LineConstraint *elems *CanRepeat
*if(CondNumEntities(int)>0)
*loop elems *OnlyInCond
*if(IsQuadratic==0)
*set var Boundary202=Operation(Boundary202+1)
mesh.boundary: *Elemsnum *cond(1) 0 0 202 *ElemsConec
*else
*set var Boundary203=Operation(Boundary203+1)
mesh.boundary: *Elemsnum *cond(1) 0 0 203 *ElemsConec
*endif
*end elems
*endif
*set var Boundary303=0
*set var Boundary306=0
*set elems(Triangle)
*set Cond SurfaceConstraint *elems *CanRepeat
*if(CondNumEntities(int)>0)
*loop elems *OnlyInCond
*if(IsQuadratic==0)
*set var Boundary303=Operation(Boundary303+1)
mesh.boundary: *Elemsnum *cond(1) 0 0 303 *ElemsConec
*else
*set var Boundary306=Operation(Boundary306+1)
mesh.boundary: *Elemsnum *cond(1) 0 0 306 *ElemsConec
*endif
*end elems
*endif
*set var Boundary404=0
*set var Boundary408=0
*set elems(Quadrilateral)
*set Cond SurfaceConstraint *elems *CanRepeat
*if(CondNumEntities(int)>0)
*loop elems *OnlyInCond
*if(IsQuadratic==0)
*set var Boundary404=Operation(Boundary404+1)
mesh.boundary: *Elemsnum *cond(1) 0 0 404 *ElemsConec
*else
*set var Boundary408=Operation(Boundary408+1)
mesh.boundary: *Elemsnum *cond(1) 0 0 408 *ElemsConec
*endif
*end elems
*endif

*set var Bulk202=0
*set var Bulk203=0
*set elems(Linear)
*loop elems
*if(ElemsMat>0)
*if(IsQuadratic==0)
*set var Bulk202=Operation(Bulk202+1)
mesh.elements: *ElemsNum *ElemsMat 202 *ElemsConec
*else
*set var Bulk203=Operation(Bulk203+1)
mesh.elements: *ElemsNum *ElemsMat 203 *ElemsConec
*endif
*endif
*end elems
*set var Bulk303=0
*set var Bulk306=0
*set elems(Triangle)
*loop elems
*if(ElemsMat>0)
*if(IsQuadratic==0)
*set var Bulk303=Operation(Bulk303+1)
mesh.elements: *Elemsnum *ElemsMat 303 *ElemsConec
*else
*set var Bulk306=Operation(Bulk306+1)
mesh.elements: *Elemsnum *ElemsMat 306 *ElemsConec
*endif
*endif
*end elems
*set var Bulk404=0
*set var Bulk408=0
*set elems(Quadrilateral)
*loop elems
*if(ElemsMat>0)
*if(IsQuadratic==0)
*set var Bulk404=Operation(Bulk404+1)
mesh.elements: *Elemsnum *ElemsMat 404 *ElemsConec
*else
*set var Bulk408=Operation(Bulk408+1)
mesh.elements: *Elemsnum *ElemsMat 408 *ElemsConec
*endif
*endif
*end elems
*set var Bulk504=0
*set var Bulk510=0
*set elems(Tetrahedra)
*loop elems
*if(ElemsMat>0)
*if(IsQuadratic==0)
*set var Bulk504=Operation(Bulk504+1)
mesh.elements: *Elemsnum *ElemsMat 504 *ElemsConec
*else
*set var Bulk510=Operation(Bulk510+1)
mesh.elements: *Elemsnum *ElemsMat 510 *ElemsConec
*endif
*endif
*end elems
*set var Bulk605=0
*set var Bulk613=0
*set elems(Pyramid)
*loop elems
*if(ElemsMat>0)
*if(IsQuadratic==0)
*set var Bulk605=Operation(Bulk605+1)
mesh.elements: *Elemsnum *ElemsMat 605 *ElemsConec
*else
*set var Bulk613=Operation(Bulk613+1)
mesh.elements: *Elemsnum *ElemsMat 613 *ElemsConec
*endif
*endif
*end elems
*set var Bulk706=0
*set var Bulk715=0
*set elems(Prism)
*loop elems
*if(ElemsMat>0)
*if(IsQuadratic==0)
*set var Bulk706=Operation(Bulk706+1)
mesh.elements: *Elemsnum *ElemsMat 706 *ElemsConec
*else
*set var Bulk715=Operation(Bulk715+1)
mesh.elements: *Elemsnum *ElemsMat 715 *ElemsConec
*endif
*endif
*end elems
*set var Bulk808=0
*set var Bulk820=0
*set elems(Hexahedra)
*loop elems
*if(ElemsMat>0)
*if(IsQuadratic==0)
*set var Bulk808=Operation(Bulk808+1)
mesh.elements: *Elemsnum *ElemsMat 808 *ElemsConec
*else
*set var Bulk820=Operation(Bulk820+1)
mesh.elements: *Elemsnum *ElemsMat 820 *ElemsConec
*endif
*endif
*end elems

*set var Boundary=Operation(Boundary202+Boundary203+Boundary303+Boundary306+Boundary404+Boundary408)
*set var Bulk=Operation(Bulk202+Bulk203+Bulk303+Bulk306+Bulk404+Bulk408+Bulk504+Bulk510+Bulk605+Bulk613+Bulk706+Bulk715+Bulk808+Bulk820)
mesh.header:*npoin *Bulk *Boundary
*set var Types=0
*if((Boundary202>0)||(Bulk202>0))
*set var Types=Operation(Types+1)
*endif
*if((Boundary203>0)||(Bulk203>0))
*set var Types=Operation(Types+1)
*endif
*if((Boundary303>0)||(Bulk303>0))
*set var Types=Operation(Types+1)
*endif
*if((Boundary306>0)||(Bulk306>0))
*set var Types=Operation(Types+1)
*endif
*if((Boundary404>0)||(Bulk404>0))
*set var Types=Operation(Types+1)
*endif
*if((Boundary408>0)||(Bulk408>0))
*set var Types=Operation(Types+1)
*endif
*if(Bulk504>0)
*set var Types=Operation(Types+1)
*endif
*if(Bulk510>0)
*set var Types=Operation(Types+1)
*endif
*if(Bulk605>0)
*set var Types=Operation(Types+1)
*endif
*if(Bulk613>0)
*set var Types=Operation(Types+1)
*endif
*if(Bulk706>0)
*set var Types=Operation(Types+1)
*endif
*if(Bulk715>0)
*set var Types=Operation(Types+1)
*endif
*if(Bulk808>0)
*set var Types=Operation(Types+1)
*endif
*if(Bulk820>0)
*set var Types=Operation(Types+1)
*endif
mesh.header:*types
*if((Boundary202>0)||(Bulk202>0))
mesh.header:    202 *Operation(Boundary202+Bulk202)
*endif
*if((Boundary203>0)||(Bulk203>0))
mesh.header:    203 *Operation(Boundary203+Bulk203)
*endif
*if((Boundary303>0)||(Bulk303>0))
mesh.header:    303 *Operation(Boundary303+Bulk303)
*endif
*if((Boundary306>0)||(Bulk306>0))
mesh.header:    306 *Operation(Boundary306+Bulk306)
*endif
*if((Boundary404>0)||(Bulk404>0))
mesh.header:    404 *Operation(Boundary404+Bulk404)
*endif
*if((Boundary408>0)||(Bulk408>0))
mesh.header:    408 *Operation(Boundary408+Bulk408)
*endif
*if(Bulk504>0)
mesh.header:    504 *Bulk504
*endif
*if(Bulk510>0)
mesh.header:    510 *Bulk510
*endif
*if(Bulk605>0)
mesh.header:    605 *Bulk605
*endif
*if(Bulk613>0)
mesh.header:    613 *Bulk613
*endif
*if(Bulk706>0)
mesh.header:    706 *Bulk706
*endif
*if(Bulk715>0)
mesh.header:    715 *Bulk715
*endif
*if(Bulk808>0)
mesh.header:    808 *Bulk808
*endif
*if(Bulk820>0)
mesh.header:    820 *Bulk820
*endif
eof

