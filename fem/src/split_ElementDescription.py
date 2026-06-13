#!/usr/bin/env python3
"""Split ElementDescription.F90 into three sub-modules + facade.

Module dependency graph (no cycles):
  ElementBasis  (no Elmer deps)
      |
  ElemInfo      (USEs ElementBasis; contains ElementInfo + EdgeElementInfo + metric)
      |
  ElementGeometry (USEs ElementBasis + ElemInfo for NormalVector→CheckPassiveElement)

Facade: USE ElementBasis + ElemInfo + ElementGeometry
"""

import os

SRC = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(SRC, 'ElementDescription.F90')

with open(ORIG, 'r') as f:
    lines = f.readlines()

LICENSE = ''.join(lines[0:43])  # lines 1-43

def get_lines(ranges):
    result = []
    for lo, hi in ranges:
        result.extend(lines[lo-1:hi])
    return result

# ---------------------------------------------------------------------------
# ElementBasis.F90
#
# Lines 83-648:   element type setup (SwapRefElemNodes, AddElementDescription,
#                 InitializeElementDescriptions, GetElementType)
# Lines 749-3036: polynomial primitives, basis dispatchers, ElementBasisDegree
# Lines 6810-6818: CrossProduct (pure vector math, no deps)
# Lines 13547-13747: GetEdgeMap (topology, no deps) + ElementDiameter
#                    (calls GetEdgeMap only)
# ---------------------------------------------------------------------------
basis_header = LICENSE + """\
#include "../config.h"

!--------------------------------------------------------------------------------
!>  Element basis functions: type setup, polynomial primitives, basis dispatchers,
!>  CrossProduct, GetEdgeMap, ElementDiameter.
!--------------------------------------------------------------------------------
MODULE ElementBasis
   USE Messages
   USE Integration
   USE LinearAlgebra
   USE CoordinateSystems
   USE PElementMaps
   USE PElementBase
   USE H1Basis
   USE Lists

   IMPLICIT NONE

   INTEGER, PARAMETER, PRIVATE :: MaxDeg  = 4, MaxDeg3 = MaxDeg**3, &
                                   MaxDeg2 = MaxDeg**2
   INTEGER, PARAMETER :: MAX_ELEMENT_NODES = 256
   LOGICAL, PRIVATE :: TypeListInitialized = .FALSE.
   TYPE(ElementType_t), POINTER :: ElementTypeList

CONTAINS

"""
basis_body = get_lines([
    (83,   583),    # SwapRefElemNodes → InitializeElementDescriptions
    (749, 3036),    # InterpolateInElement1D → ElementBasisDegree
    (6810, 6818),   # CrossProduct
    (13301, 13359), # InterpolateInElement general dispatcher (only calls 1D/2D/3D variants)
    (13547, 13747), # GetEdgeMap + ElementDiameter
])
basis_footer = "\nEND MODULE ElementBasis\n"

with open(os.path.join(SRC, 'ElementBasis.F90'), 'w') as f:
    f.write(basis_header)
    f.writelines(basis_body)
    f.write(basis_footer)

print("ElementBasis.F90 written:", sum(1 for _ in basis_body), "body lines")

# ---------------------------------------------------------------------------
# ElemInfo.F90  (MODULE name avoids collision with FUNCTION ElementInfo inside)
#
# Lines 649-748:    StabParam (calls ElementInfo + ElementDiameter from ElementBasis)
# Lines 3040-6809:  EdgeElementStyle, ElementInfo, ElementInfoVec family,
#                   ElementSize, FaceElementInfo, PiolaTransformationData,
#                   FaceElementOrientation, FaceElementBasisOrdering, PickActiveFace
# Lines 6819-13541: separators + EdgeElementInfo, DOFs helpers, WeightedWhitneyForms,
#                   GetEdgeBasis, mGetElementDOFs, CheckMetric (#ifdef HAVE_QP),
#                   ElementMetric, ElementMetricQP, ElementMetricVec,
#                   GlobalFirstDerivativesInternal, GlobalFirstDerivatives,
#                   InterpolateInElement, GlobalSecondDerivatives
# ---------------------------------------------------------------------------
info_header = LICENSE + """\
#include "../config.h"

!--------------------------------------------------------------------------------
!>  Element information assembly: ElementInfo, EdgeElementInfo, FaceElementInfo,
!>  metric computation, global derivative computation. StabParam also lives here.
!--------------------------------------------------------------------------------
MODULE ElemInfo
   USE ElementBasis
   USE ElementGeometry
   USE Messages
   USE Integration
   USE LinearAlgebra
   USE CoordinateSystems
   USE PElementMaps
   USE PElementBase
   USE H1Basis
   USE Lists

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: GetElementType, StabParam, &
             EdgeElementStyle, ElementInfo, ElementInfoVec, &
             ElementInfoVec_ComputePElementBasis, ElementInfoVec_ElementBasisToGlobal, &
             ElementSize, &
             FaceElementInfo, PiolaTransformationData, &
             FaceElementOrientation, FaceElementBasisOrdering, PickActiveFace, &
             EdgeElementInfo, &
             TriangleFaceDofsOrdering, TriangleFaceDofsOrdering2nd, &
             TriangleFaceDofsOrdering2, SquareFaceDofsOrdering, &
             ReorderingAndSignReversionsData, &
             EdgeWhitneyComponents2D, FaceWhitneyComponents2D, WeightedWhitneyForms, &
             GetEdgeBasis, mGetElementDOFs, &
             CheckMetric, ElementMetric, ElementMetricQP, ElementMetricVec, &
             GlobalFirstDerivativesInternal, GlobalFirstDerivatives, &
             GlobalSecondDerivatives

CONTAINS

"""
info_body = get_lines([
    (584,   748),   # GetElementType (calls StabParam) + StabParam
    (3040,  6809),  # EdgeElementStyle → PickActiveFace
    (6819, 13300),  # EdgeElementInfo → InterpolateInElement (excluded, moved to ElementBasis)
    (13360, 13541), # GlobalSecondDerivatives (after InterpolateInElement)
])
info_footer = "\nEND MODULE ElemInfo\n"

with open(os.path.join(SRC, 'ElemInfo.F90'), 'w') as f:
    f.write(info_header)
    f.writelines(info_body)
    f.write(info_footer)

print("ElemInfo.F90 written:", sum(1 for _ in info_body), "body lines")

# ---------------------------------------------------------------------------
# ElementGeometry.F90
#
# Lines 13748-15651: TriangleInside, QuadInside, TetraInside, BrickInside,
#                    CheckPassiveElement, CheckNormalDirection,
#                    CheckNormalDirectionParent, NormalVector, NormalVectorLinear,
#                    SurfaceVector, LineFaceIntersection, LineFaceIntersection2,
#                    PointFaceDistance, GlobalToLocal, getTriangleFaceDirection,
#                    getSquareFaceDirection, wedgeOrdering, ComputeRotationMatrix,
#                    CutSingleElement, SplitSingleElement
# ---------------------------------------------------------------------------
geom_header = LICENSE + """\
#include "../config.h"

!--------------------------------------------------------------------------------
!>  Element geometry: normal/surface vectors, coordinate mapping, inside tests,
!>  cutting/splitting.
!--------------------------------------------------------------------------------
MODULE ElementGeometry
   USE ElementBasis
   USE Messages
   USE Integration
   USE LinearAlgebra
   USE CoordinateSystems
   USE PElementMaps
   USE PElementBase
   USE Lists

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: TriangleInside, QuadInside, TetraInside, BrickInside, &
             CheckPassiveElement, CheckNormalDirection, CheckNormalDirectionParent, &
             NormalVector, SurfaceVector, &
             LineFaceIntersection, LineFaceIntersection2, PointFaceDistance, &
             GlobalToLocal, &
             getTriangleFaceDirection, getSquareFaceDirection, wedgeOrdering, &
             ComputeRotationMatrix, CutSingleElement, SplitSingleElement

CONTAINS

"""
geom_body = get_lines([(13748, 15651)])
geom_footer = "\nEND MODULE ElementGeometry\n"

with open(os.path.join(SRC, 'ElementGeometry.F90'), 'w') as f:
    f.write(geom_header)
    f.writelines(geom_body)
    f.write(geom_footer)

print("ElementGeometry.F90 written:", sum(1 for _ in geom_body), "body lines")

# ---------------------------------------------------------------------------
# ElementDescription.F90 — facade re-exporting all sub-modules
# ---------------------------------------------------------------------------
facade = LICENSE + """\
#include "../config.h"

!--------------------------------------------------------------------------------
!>  Facade module: re-exports ElementBasis, ElemInfo, and ElementGeometry so
!>  that existing 'USE ElementDescription' continues to work unchanged.
!--------------------------------------------------------------------------------
!> \\ingroup ElmerLib
!> \\{
MODULE ElementDescription
   USE ElementBasis
   USE ElemInfo
   USE ElementGeometry
   IMPLICIT NONE
END MODULE ElementDescription
!> \\}
"""

with open(os.path.join(SRC, 'ElementDescription.F90'), 'w') as f:
    f.write(facade)

print("ElementDescription.F90 (facade) written")
print("Done.")
