// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_samrai_interface
#define included_IBTK_samrai_interface

#include <tbox/Pointer.h>

namespace IBTK
{
template <class T>
using SAMRAIPointer = SAMRAI::tbox::Pointer<T>;

template <class T, class... Args>
inline SAMRAIPointer<T>
make_samrai_shared(Args&&... args)
{
    return SAMRAIPointer<T>(new T(std::forward<Args>(args)...));
}

template <class T>
inline SAMRAIPointer<T>
make_samrai_shared(std::size_t N)
{
    return SAMRAIPointer<T>(new T[N]);
}

} // namespace IBTK

namespace SAMRAI
{
namespace algs
{
template <int DIM>
class HyperbolicLevelIntegrator;
using HyperbolicLevelIntegratorNd = HyperbolicLevelIntegrator<NDIM>;

template <int DIM>
class HyperbolicPatchStrategy;
using HyperbolicPatchStrategyNd = HyperbolicPatchStrategy<NDIM>;

template <int DIM>
class TimeRefinementIntegrator;
using TimeRefinementIntegratorNd = TimeRefinementIntegrator<NDIM>;
} // namespace algs

namespace appu
{
template <int DIM>
class VisItDataWriter;
using VisItDataWriterNd = VisItDataWriter<NDIM>;
} // namespace appu

namespace geom
{
template <int DIM>
class CartesianCellDoubleLinearRefine;
using CartesianCellDoubleLinearRefineNd = CartesianCellDoubleLinearRefine<NDIM>;

template <int DIM>
class CartesianCellDoubleConservativeLinearRefine;
using CartesianCellDoubleConservativeLinearRefineNd = CartesianCellDoubleConservativeLinearRefine<NDIM>;

template <int DIM>
class CartesianCellDoubleWeightedAverage;
using CartesianCellDoubleWeightedAverageNd = CartesianCellDoubleWeightedAverage<NDIM>;

template <int DIM>
class CartesianGridGeometry;
using CartesianGridGeometryNd = CartesianGridGeometry<NDIM>;

template <int DIM>
class CartesianPatchGeometry;
using CartesianPatchGeometryNd = CartesianPatchGeometry<NDIM>;

template <int DIM>
class CartesianSideDoubleConservativeLinearRefine;
using CartesianSideDoubleConservativeLinearRefineNd = CartesianSideDoubleConservativeLinearRefine<NDIM>;

template <int DIM>
class CartesianSideDoubleWeightedAverage;
using CartesianSideDoubleWeightedAverageNd = CartesianSideDoubleWeightedAverage<NDIM>;
} // namespace geom

namespace hier
{
template <int DIM>
class BasePatchHierarchy;
using BasePatchHierarchyNd = BasePatchHierarchy<NDIM>;

template <int DIM>
class BasePatchLevel;
using BasePatchLevelNd = BasePatchLevel<NDIM>;

template <int DIM>
class BoundaryBox;
using BoundaryBoxNd = BoundaryBox<NDIM>;

template <int DIM>
class Box;
using BoxNd = Box<NDIM>;

template <int DIM>
class BoxArray;
using BoxArrayNd = BoxArray<NDIM>;

template <int DIM>
class BoxGeometry;
using BoxGeometryNd = BoxGeometry<NDIM>;

template <int DIM>
class BoxIterator;
using BoxIteratorNd = BoxIterator<NDIM>;

template <int DIM>
class BoxList;
using BoxListNd = BoxList<NDIM>;

template <int DIM>
class BoxOverlap;
using BoxOverlapNd = BoxOverlap<NDIM>;

template <int DIM>
class BoxTree;
using BoxTreeNd = BoxTree<NDIM>;

template <int DIM>
class BoundaryLookupTable;
using BoundaryLookupTableNd = BoundaryLookupTable<NDIM>;

template <int DIM>
class CoarseFineBoundary;
using CoarseFineBoundaryNd = CoarseFineBoundary<NDIM>;

template <int DIM>
class GridGeometry;
using GridGeometryNd = GridGeometry<NDIM>;

template <int DIM>
class Index;
using IndexNd = Index<NDIM>;

template <int DIM>
class IntVector;
using IntVectorNd = IntVector<NDIM>;

template <int DIM>
class Patch;
using PatchNd = Patch<NDIM>;

template <int DIM>
class PatchData;
using PatchDataNd = PatchData<NDIM>;

template <int DIM>
class PatchDataFactory;
using PatchDataFactoryNd = PatchDataFactory<NDIM>;

template <int DIM>
class PatchDescriptor;
using PatchDescriptorNd = PatchDescriptor<NDIM>;

template <int DIM>
class PatchGeometry;
using PatchGeometryNd = PatchGeometry<NDIM>;

template <int DIM>
class PatchHierarchy;
using PatchHierarchyNd = PatchHierarchy<NDIM>;

template <int DIM>
class PatchLevel;
using PatchLevelNd = PatchLevel<NDIM>;

template <int DIM>
class Variable;
using VariableNd = Variable<NDIM>;

template <int DIM>
class VariableDatabase;
using VariableDatabaseNd = VariableDatabase<NDIM>;
} // namespace hier

namespace math
{
template <int DIM, class T>
class ArrayDataBasicOps;
template <class T>
using ArrayDataBasicOpsNd = ArrayDataBasicOps<NDIM, T>;

template <int DIM, class T>
class ArrayDataMiscellaneousOpsReal;
template <class T>
using ArrayDataMiscellaneousOpsRealNd = ArrayDataMiscellaneousOpsReal<NDIM, T>;

template <int DIM, class T>
class ArrayDataNormOpsComplex;
template <class T>
using ArrayDataNormOpsComplexNd = ArrayDataNormOpsComplex<NDIM, T>;

template <int DIM>
class ArrayDataNormOpsInteger;
using ArrayDataNormOpsIntegerNd = ArrayDataNormOpsInteger<NDIM>;

template <int DIM, class T>
class ArrayDataNormOpsReal;
template <class T>
using ArrayDataNormOpsRealNd = ArrayDataNormOpsReal<NDIM, T>;

template <int DIM>
class HierarchyDataOpsManager;
using HierarchyDataOpsManagerNd = HierarchyDataOpsManager<NDIM>;

template <int DIM>
class HierarchyCellDataOpsInteger;
using HierarchyCellDataOpsIntegerNd = HierarchyCellDataOpsInteger<NDIM>;

template <int DIM, class T>
class HierarchyCellDataOpsReal;
template <class T>
using HierarchyCellDataOpsRealNd = HierarchyCellDataOpsReal<NDIM, T>;

template <int DIM>
class HierarchyDataOpsInteger;
using HierarchyDataOpsIntegerNd = HierarchyDataOpsInteger<NDIM>;

template <int DIM, class T>
class HierarchyDataOpsReal;
template <class T>
using HierarchyDataOpsRealNd = HierarchyDataOpsReal<NDIM, T>;

template <int DIM>
class HierarchyEdgeDataOpsInteger;
using HierarchyEdgeDataOpsIntegerNd = HierarchyEdgeDataOpsInteger<NDIM>;

template <int DIM, class T>
class HierarchyEdgeDataOpsReal;
template <class T>
using HierarchyEdgeDataOpsRealNd = HierarchyEdgeDataOpsReal<NDIM, T>;

template <int DIM>
class HierarchyFaceDataOpsInteger;
using HierarchyFaceDataOpsIntegerNd = HierarchyFaceDataOpsInteger<NDIM>;

template <int DIM, class T>
class HierarchyFaceDataOpsReal;
template <class T>
using HierarchyFaceDataOpsRealNd = HierarchyFaceDataOpsReal<NDIM, T>;

template <int DIM>
class HierarchyNodeDataOpsInteger;
using HierarchyNodeDataOpsIntegerNd = HierarchyNodeDataOpsInteger<NDIM>;

template <int DIM, class T>
class HierarchyNodeDataOpsReal;
template <class T>
using HierarchyNodeDataOpsRealNd = HierarchyNodeDataOpsReal<NDIM, T>;

template <int DIM>
class HierarchySideDataOpsInteger;
using HierarchySideDataOpsIntegerNd = HierarchySideDataOpsInteger<NDIM>;

template <int DIM, class T>
class HierarchySideDataOpsReal;
template <class T>
using HierarchySideDataOpsRealNd = HierarchySideDataOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchCellDataBasicOps;
template <class T>
using PatchCellDataBasicOpsNd = PatchCellDataBasicOps<NDIM, T>;

template <int DIM, class T>
class PatchCellDataMiscellaneousOpsReal;
template <class T>
using PatchCellDataMiscellaneousOpsRealNd = PatchCellDataMiscellaneousOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchCellDataNormOpsComplex;
template <class T>
using PatchCellDataNormOpsComplexNd = PatchCellDataNormOpsComplex<NDIM, T>;

template <int DIM, class T>
class PatchCellDataNormOpsReal;
template <class T>
using PatchCellDataNormOpsRealNd = PatchCellDataNormOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchCellDataOpsComplex;
template <class T>
using PatchCellDataOpsComplexNd = PatchCellDataOpsComplex<NDIM, T>;

template <int DIM>
class PatchCellDataOpsInteger;
using PatchCellDataOpsIntegerNd = PatchCellDataOpsInteger<NDIM>;

template <int DIM, class T>
class PatchCellDataOpsReal;
template <class T>
using PatchCellDataOpsRealNd = PatchCellDataOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchEdgeDataBasicOps;
template <class T>
using PatchEdgeDataBasicOpsNd = PatchEdgeDataBasicOps<NDIM, T>;

template <int DIM, class T>
class PatchEdgeDataMisEdgeaneousOpsReal;
template <class T>
using PatchEdgeDataMisEdgeaneousOpsRealNd = PatchEdgeDataMisEdgeaneousOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchEdgeDataNormOpsComplex;
template <class T>
using PatchEdgeDataNormOpsComplexNd = PatchEdgeDataNormOpsComplex<NDIM, T>;

template <int DIM, class T>
class PatchEdgeDataNormOpsReal;
template <class T>
using PatchEdgeDataNormOpsRealNd = PatchEdgeDataNormOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchEdgeDataOpsComplex;
template <class T>
using PatchEdgeDataOpsComplexNd = PatchEdgeDataOpsComplex<NDIM, T>;

template <int DIM>
class PatchEdgeDataOpsInteger;
using PatchEdgeDataOpsIntegerNd = PatchEdgeDataOpsInteger<NDIM>;

template <int DIM, class T>
class PatchEdgeDataOpsReal;
template <class T>
using PatchEdgeDataOpsRealNd = PatchEdgeDataOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchFaceDataBasicOps;
template <class T>
using PatchFaceDataBasicOpsNd = PatchFaceDataBasicOps<NDIM, T>;

template <int DIM, class T>
class PatchFaceDataMisFaceaneousOpsReal;
template <class T>
using PatchFaceDataMisFaceaneousOpsRealNd = PatchFaceDataMisFaceaneousOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchFaceDataNormOpsComplex;
template <class T>
using PatchFaceDataNormOpsComplexNd = PatchFaceDataNormOpsComplex<NDIM, T>;

template <int DIM, class T>
class PatchFaceDataNormOpsReal;
template <class T>
using PatchFaceDataNormOpsRealNd = PatchFaceDataNormOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchFaceDataOpsComplex;
template <class T>
using PatchFaceDataOpsComplexNd = PatchFaceDataOpsComplex<NDIM, T>;

template <int DIM>
class PatchFaceDataOpsInteger;
using PatchFaceDataOpsIntegerNd = PatchFaceDataOpsInteger<NDIM>;

template <int DIM, class T>
class PatchFaceDataOpsReal;
template <class T>
using PatchFaceDataOpsRealNd = PatchFaceDataOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchNodeDataBasicOps;
template <class T>
using PatchNodeDataBasicOpsNd = PatchNodeDataBasicOps<NDIM, T>;

template <int DIM, class T>
class PatchNodeDataMisNodeaneousOpsReal;
template <class T>
using PatchNodeDataMisNodeaneousOpsRealNd = PatchNodeDataMisNodeaneousOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchNodeDataNormOpsComplex;
template <class T>
using PatchNodeDataNormOpsComplexNd = PatchNodeDataNormOpsComplex<NDIM, T>;

template <int DIM, class T>
class PatchNodeDataNormOpsReal;
template <class T>
using PatchNodeDataNormOpsRealNd = PatchNodeDataNormOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchNodeDataOpsComplex;
template <class T>
using PatchNodeDataOpsComplexNd = PatchNodeDataOpsComplex<NDIM, T>;

template <int DIM>
class PatchNodeDataOpsInteger;
using PatchNodeDataOpsIntegerNd = PatchNodeDataOpsInteger<NDIM>;

template <int DIM, class T>
class PatchNodeDataOpsReal;
template <class T>
using PatchNodeDataOpsRealNd = PatchNodeDataOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchSideDataBasicOps;
template <class T>
using PatchSideDataBasicOpsNd = PatchSideDataBasicOps<NDIM, T>;

template <int DIM, class T>
class PatchSideDataMisSideaneousOpsReal;
template <class T>
using PatchSideDataMisSideaneousOpsRealNd = PatchSideDataMisSideaneousOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchSideDataNormOpsComplex;
template <class T>
using PatchSideDataNormOpsComplexNd = PatchSideDataNormOpsComplex<NDIM, T>;

template <int DIM, class T>
class PatchSideDataNormOpsReal;
template <class T>
using PatchSideDataNormOpsRealNd = PatchSideDataNormOpsReal<NDIM, T>;

template <int DIM, class T>
class PatchSideDataOpsComplex;
template <class T>
using PatchSideDataOpsComplexNd = PatchSideDataOpsComplex<NDIM, T>;

template <int DIM>
class PatchSideDataOpsInteger;
using PatchSideDataOpsIntegerNd = PatchSideDataOpsInteger<NDIM>;

template <int DIM, class T>
class PatchSideDataOpsReal;
template <class T>
using PatchSideDataOpsRealNd = PatchSideDataOpsReal<NDIM, T>;
} // namespace math

namespace mesh
{
template <int DIM>
class BergerRigoutsos;
using BergerRigoutsosNd = BergerRigoutsos<NDIM>;

template <int DIM>
class BoxGeneratorStrategy;
using BoxGeneratorStrategyNd = BoxGeneratorStrategy<NDIM>;

template <int DIM>
class GriddingAlgorithm;
using GriddingAlgorithmNd = GriddingAlgorithm<NDIM>;

template <int DIM>
class LoadBalancer;
using LoadBalancerNd = LoadBalancer<NDIM>;

template <int DIM>
class StandardTagAndInitStrategy;
using StandardTagAndInitStrategyNd = StandardTagAndInitStrategy<NDIM>;

template <int DIM>
class StandardTagAndInitialize;
using StandardTagAndInitializeNd = StandardTagAndInitialize<NDIM>;

template <int DIM>
class TagAndInitializeStrategy;
using TagAndInitializeStrategyNd = TagAndInitializeStrategy<NDIM>;
} // namespace mesh

namespace pdat
{
template <int DIM, class T>
class ArrayData;
template <class T>
using ArrayDataNd = ArrayData<NDIM, T>;

template <int DIM, class T>
class CellData;
template <class T>
using CellDataNd = CellData<NDIM, T>;

template <int DIM, class T>
class CellDataFactory;
template <class T>
using CellDataFactoryNd = CellDataFactory<NDIM, T>;

template <int DIM>
class CellDoubleConstantRefine;
using CellDoubleConstantRefineNd = CellDoubleConstantRefine<NDIM>;

template <int DIM>
class CellGeometry;
using CellGeometryNd = CellGeometry<NDIM>;

template <int DIM>
class CellIndex;
using CellIndexNd = CellIndex<NDIM>;

template <int DIM>
class CellIterator;
using CellIteratorNd = CellIterator<NDIM>;

template <int DIM>
class CellOverlap;
using CellOverlapNd = CellOverlap<NDIM>;

template <int DIM, class T>
class CellVariable;
template <class T>
using CellVariableNd = CellVariable<NDIM, T>;

template <int DIM, class T>
class EdgeData;
template <class T>
using EdgeDataNd = EdgeData<NDIM, T>;

template <int DIM, class T>
class EdgeDataFactory;
template <class T>
using EdgeDataFactoryNd = EdgeDataFactory<NDIM, T>;

template <int DIM>
class EdgeGeometry;
using EdgeGeometryNd = EdgeGeometry<NDIM>;

template <int DIM>
class EdgeIndex;
using EdgeIndexNd = EdgeIndex<NDIM>;

template <int DIM>
class EdgeIterator;
using EdgeIteratorNd = EdgeIterator<NDIM>;

template <int DIM>
class EdgeOverlap;
using EdgeOverlapNd = EdgeOverlap<NDIM>;

template <int DIM, class T>
class EdgeVariable;
template <class T>
using EdgeVariableNd = EdgeVariable<NDIM, T>;

template <int DIM, class T>
class FaceData;
template <class T>
using FaceDataNd = FaceData<NDIM, T>;

template <int DIM, class T>
class FaceDataFactory;
template <class T>
using FaceDataFactoryNd = FaceDataFactory<NDIM, T>;

template <int DIM>
class FaceGeometry;
using FaceGeometryNd = FaceGeometry<NDIM>;

template <int DIM>
class FaceIndex;
using FaceIndexNd = FaceIndex<NDIM>;

template <int DIM>
class FaceIterator;
using FaceIteratorNd = FaceIterator<NDIM>;

template <int DIM>
class FaceOverlap;
using FaceOverlapNd = FaceOverlap<NDIM>;

template <int DIM, class T>
class FaceVariable;
template <class T>
using FaceVariableNd = FaceVariable<NDIM, T>;

template <int DIM, class T, class BOX_GEOMETRY>
class IndexData;
template <class T, class BOX_GEOMETRY>
using IndexDataNd = IndexData<NDIM, T, BOX_GEOMETRY>;

template <int DIM, class T, class BOX_GEOMETRY>
class IndexDataFactory;
template <class T, class BOX_GEOMETRY>
using IndexDataFactoryNd = IndexDataFactory<NDIM, T, BOX_GEOMETRY>;

template <int DIM, class T, class BOX_GEOMETRY>
class IndexDataNode;
template <class T, class BOX_GEOMETRY>
using IndexDataNodeNd = IndexDataNode<NDIM, T, BOX_GEOMETRY>;

template <int DIM, class T, class BOX_GEOMETRY>
class IndexIterator;
template <class T, class BOX_GEOMETRY>
using IndexIteratorNd = IndexIterator<NDIM, T, BOX_GEOMETRY>;

template <int DIM, class T, class BOX_GEOMETRY>
class IndexVariable;
template <class T, class BOX_GEOMETRY>
using IndexVariableNd = IndexVariable<NDIM, T, BOX_GEOMETRY>;

template <int DIM, class T>
class NodeData;
template <class T>
using NodeDataNd = NodeData<NDIM, T>;

template <int DIM, class T>
class NodeDataFactory;
template <class T>
using NodeDataFactoryNd = NodeDataFactory<NDIM, T>;

template <int DIM>
class NodeGeometry;
using NodeGeometryNd = NodeGeometry<NDIM>;

template <int DIM>
class NodeIndex;
using NodeIndexNd = NodeIndex<NDIM>;

template <int DIM>
class NodeIterator;
using NodeIteratorNd = NodeIterator<NDIM>;

template <int DIM>
class NodeOverlap;
using NodeOverlapNd = NodeOverlap<NDIM>;

template <int DIM, class T>
class NodeVariable;
template <class T>
using NodeVariableNd = NodeVariable<NDIM, T>;

template <int DIM, class T>
class OuteredgeData;
template <class T>
using OuteredgeDataNd = OuteredgeData<NDIM, T>;

template <int DIM, class T>
class OuteredgeDataFactory;
template <class T>
using OuteredgeDataFactoryNd = OuteredgeDataFactory<NDIM, T>;

template <int DIM>
class OuteredgeGeometry;
using OuteredgeGeometryNd = OuteredgeGeometry<NDIM>;

template <int DIM>
class OuteredgeIndex;
using OuteredgeIndexNd = OuteredgeIndex<NDIM>;

template <int DIM>
class OuteredgeIterator;
using OuteredgeIteratorNd = OuteredgeIterator<NDIM>;

template <int DIM, class T>
class OuteredgeVariable;
template <class T>
using OuteredgeVariableNd = OuteredgeVariable<NDIM, T>;

template <int DIM, class T>
class OuterfaceData;
template <class T>
using OuterfaceDataNd = OuterfaceData<NDIM, T>;

template <int DIM, class T>
class OuterfaceDataFactory;
template <class T>
using OuterfaceDataFactoryNd = OuterfaceDataFactory<NDIM, T>;

template <int DIM>
class OuterfaceGeometry;
using OuterfaceGeometryNd = OuterfaceGeometry<NDIM>;

template <int DIM>
class OuterfaceIndex;
using OuterfaceIndexNd = OuterfaceIndex<NDIM>;

template <int DIM>
class OuterfaceIterator;
using OuterfaceIteratorNd = OuterfaceIterator<NDIM>;

template <int DIM, class T>
class OuterfaceVariable;
template <class T>
using OuterfaceVariableNd = OuterfaceVariable<NDIM, T>;

template <int DIM, class T>
class OuternodeData;
template <class T>
using OuternodeDataNd = OuternodeData<NDIM, T>;

template <int DIM, class T>
class OuternodeDataFactory;
template <class T>
using OuternodeDataFactoryNd = OuternodeDataFactory<NDIM, T>;

template <int DIM>
class OuternodeGeometry;
using OuternodeGeometryNd = OuternodeGeometry<NDIM>;

template <int DIM>
class OuternodeIndex;
using OuternodeIndexNd = OuternodeIndex<NDIM>;

template <int DIM>
class OuternodeIterator;
using OuternodeIteratorNd = OuternodeIterator<NDIM>;

template <int DIM, class T>
class OuternodeVariable;
template <class T>
using OuternodeVariableNd = OuternodeVariable<NDIM, T>;

template <int DIM, class T>
class OutersideData;
template <class T>
using OutersideDataNd = OutersideData<NDIM, T>;

template <int DIM, class T>
class OutersideDataFactory;
template <class T>
using OutersideDataFactoryNd = OutersideDataFactory<NDIM, T>;

template <int DIM>
class OutersideGeometry;
using OutersideGeometryNd = OutersideGeometry<NDIM>;

template <int DIM>
class OutersideIndex;
using OutersideIndexNd = OutersideIndex<NDIM>;

template <int DIM>
class OutersideIterator;
using OutersideIteratorNd = OutersideIterator<NDIM>;

template <int DIM, class T>
class OutersideVariable;
template <class T>
using OutersideVariableNd = OutersideVariable<NDIM, T>;

template <int DIM, class T>
class SideData;
template <class T>
using SideDataNd = SideData<NDIM, T>;

template <int DIM, class T>
class SideDataFactory;
template <class T>
using SideDataFactoryNd = SideDataFactory<NDIM, T>;

template <int DIM>
class SideGeometry;
using SideGeometryNd = SideGeometry<NDIM>;

template <int DIM>
class SideIndex;
using SideIndexNd = SideIndex<NDIM>;

template <int DIM>
class SideIterator;
using SideIteratorNd = SideIterator<NDIM>;

template <int DIM>
class SideOverlap;
using SideOverlapNd = SideOverlap<NDIM>;

template <int DIM, class T>
class SideVariable;
template <class T>
using SideVariableNd = SideVariable<NDIM, T>;
} // namespace pdat

namespace solv
{
template <int DIM>
class LocationIndexRobinBcCoefs;
using LocationIndexRobinBcCoefsNd = LocationIndexRobinBcCoefs<NDIM>;

template <int DIM>
class RobinBcCoefStrategy;
using RobinBcCoefStrategyNd = RobinBcCoefStrategy<NDIM>;

template <int DIM, class T>
class SAMRAIVectorReal;
template <class T>
using SAMRAIVectorRealNd = SAMRAIVectorReal<NDIM, T>;
} // namespace solv

namespace xfer
{
template <int DIM>
class BoxGeometryFillPattern;
using BoxGeometryFillPatternNd = BoxGeometryFillPattern<NDIM>;

template <int DIM>
class CoarsenAlgorithm;
using CoarsenAlgorithmNd = CoarsenAlgorithm<NDIM>;

template <int DIM>
class CoarsenOperator;
using CoarsenOperatorNd = CoarsenOperator<NDIM>;

template <int DIM>
class CoarsenPatchStrategy;
using CoarsenPatchStrategyNd = CoarsenPatchStrategy<NDIM>;

template <int DIM>
class CoarsenSchedule;
using CoarsenScheduleNd = CoarsenSchedule<NDIM>;

template <int DIM>
class Geometry;
using GeometryNd = Geometry<NDIM>;

template <int DIM>
class RefineAlgorithm;
using RefineAlgorithmNd = RefineAlgorithm<NDIM>;

template <int DIM>
class RefineClasses;
using RefineClassesNd = RefineClasses<NDIM>;

template <int DIM>
class RefineOperator;
using RefineOperatorNd = RefineOperator<NDIM>;

template <int DIM>
class RefinePatchStrategy;
using RefinePatchStrategyNd = RefinePatchStrategy<NDIM>;

template <int DIM>
class RefineSchedule;
using RefineScheduleNd = RefineSchedule<NDIM>;

template <int DIM>
class VariableFillPattern;
using VariableFillPatternNd = VariableFillPattern<NDIM>;
} // namespace xfer
} // namespace SAMRAI

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_samrai_interface
