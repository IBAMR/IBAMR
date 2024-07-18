# algs
echo "algs"
perl -pi -e 's/HyperbolicLevelIntegrator\<NDIM\>/HyperbolicLevelIntegratorNd/g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/HyperbolicPatchStrategy\<NDIM\>/HyperbolicPatchStrategyNd/g'     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/TimeRefinementIntegrator\<NDIM\>/TimeRefinementIntegratorNd/g'   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

# appu
echo "appu"
perl -pi -e 's/VisItDataWriter\<NDIM\>/VisItDataWriterNd/g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

# geom
echo "geom"
perl -pi -e 's/CartesianCellDoubleConservativeLinearRefine\<NDIM\>/CartesianCellDoubleConservativeLinearRefineNd/g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CartesianCellDoubleLinearRefine\<NDIM\>/CartesianCellDoubleLinearRefineNd/g'                         `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CartesianCellDoubleWeightedAverage\<NDIM\>/CartesianCellDoubleWeightedAverageNd/g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CartesianGridGeometry\<NDIM\>/CartesianGridGeometryNd/g'                                             `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CartesianPatchGeometry\<NDIM\>/CartesianPatchGeometryNd/g'                                           `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CartesianSideDoubleConservativeLinearRefine\<NDIM\>/CartesianSideDoubleConservativeLinearRefineNd/g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CartesianSideDoubleWeightedAverage\<NDIM\>/CartesianSideDoubleWeightedAverageNd/g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

# hier
echo "hier"
perl -pi -e 's/Box\<NDIM\>/BoxNd/g'                                 `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/BoxArray\<NDIM\>/BoxArrayNd/g'                       `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/BoxList\<NDIM\>/BoxListNd/g'                         `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/BoxGeometry\<NDIM\>/BoxGeometryNd/g'                 `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/BoxIterator\<NDIM\>/BoxIteratorNd/g'                 `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/BoxOverlap\<NDIM\>/BoxOverlapNd/g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/BoxTree\<NDIM\>/BoxTreeNd/g'                         `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/BoxTree\<NDIM\>/BoxTreeNd/g'                         `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/BoundaryLookupTable\<NDIM\>/BoundaryLookupTableNd/g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/GridGeometry\<NDIM\>/GridGeometryNd/g'               `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/Index\<NDIM\>/IndexNd/g'                             `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/IntVector\<NDIM\>/IntVectorNd/g'                     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/Patch\<NDIM\>/PatchNd/g'                             `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchData\<NDIM\>/PatchDataNd/g'                     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchDataFactory\<NDIM\>/PatchDataFactoryNd/g'       `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchDescriptor\<NDIM\>/PatchDescriptorNd/g'         `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchGeometry\<NDIM\>/PatchGeometryNd/g'             `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchHierarchy\<NDIM\>/PatchHierarchyNd/g'           `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchLevel\<NDIM\>/PatchLevelNd/g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/Variable\<NDIM\>/VariableNd/g'                       `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/VariableDatabase\<NDIM\>/VariableDatabaseNd/g'       `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

# math
echo "math"

echo "  ArrayData"
perl -pi -e 's/ArrayDataBasicOps\<NDIM\,/ArrayDataBasicOpsNd\</g'                                 `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/ArrayDataMiscellaneousOpsReal\<NDIM\,/ArrayDataMiscellaneousOpsRealNd\</g'         `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/ArrayDataNormOpsComplex\<NDIM\,/ArrayDataNormOpsComplexNd\</g'                     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/ArrayDataNormOpsInteger\<NDIM\,/ArrayDataNormOpsIntegerNd\</g'                     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/ArrayDataNormOpsReal\<NDIM\,/ArrayDataNormOpsRealNd\</g'                           `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

echo "  HierarchyData"
perl -pi -e 's/HierarchyCellDataOpsInteger\<NDIM\>/HierarchyCellDataOpsIntegerNd/g'               `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/HierarchyCellDataOpsReal\<NDIM\,/HierarchyCellDataOpsRealNd\</g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/HierarchyDataOpsManager\<NDIM\>/HierarchyDataOpsManagerNd/g'                       `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/HierarchyDataOpsInteger\<NDIM\>/HierarchyDataOpsIntegerNd/g'                       `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/HierarchyDataOpsReal\<NDIM\,/HierarchyDataOpsRealNd\</g'                           `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/HierarchyEdgeDataOpsInteger\<NDIM\>/HierarchyEdgeDataOpsIntegerNd/g'               `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/HierarchyEdgeDataOpsReal\<NDIM\,/HierarchyEdgeDataOpsRealNd\</g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/HierarchyFaceDataOpsInteger\<NDIM\>/HierarchyFaceDataOpsIntegerNd/g'               `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/HierarchyFaceDataOpsReal\<NDIM\,/HierarchyFaceDataOpsRealNd\</g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/HierarchyNodeDataOpsInteger\<NDIM\>/HierarchyNodeDataOpsIntegerNd/g'               `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/HierarchyNodeDataOpsReal\<NDIM\,/HierarchyNodeDataOpsRealNd\</g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/HierarchySideDataOpsInteger\<NDIM\>/HierarchySideDataOpsIntegerNd/g'               `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/HierarchySideDataOpsReal\<NDIM\,/HierarchySideDataOpsRealNd\</g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

echo "  PatchCellData"
perl -pi -e 's/PatchCellDataBasicOps\<NDIM\,/PatchCellDataBasicOpsNd\</g'                         `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchCellDataMiscellaneousOpsReal\<NDIM\,/PatchCellDataMiscellaneousOpsRealNd\</g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchCellDataNormOpsComplex\<NDIM\,/PatchCellDataNormOpsComplexNd\</g'             `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchCellDataNormOpsReal\<NDIM\,/PatchCellDataNormOpsRealNd\</g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchCellDataOpsComplex\<NDIM\,/PatchCellDataOpsComplexNd\</g'                     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchCellDataOpsInteger\<NDIM\,/PatchCellDataOpsIntegerNd\</g'                     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchCellDataOpsReal\<NDIM\,/PatchCellDataOpsRealNd\</g'                           `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

echo "  PatchEdgeData"
perl -pi -e 's/PatchEdgeDataBasicOps\<NDIM\,/PatchEdgeDataBasicOpsNd\</g'                         `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchEdgeDataMisEdgeaneousOpsReal\<NDIM\,/PatchEdgeDataMisEdgeaneousOpsRealNd\</g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchEdgeDataNormOpsComplex\<NDIM\,/PatchEdgeDataNormOpsComplexNd\</g'             `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchEdgeDataNormOpsReal\<NDIM\,/PatchEdgeDataNormOpsRealNd\</g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchEdgeDataOpsComplex\<NDIM\,/PatchEdgeDataOpsComplexNd\</g'                     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchEdgeDataOpsInteger\<NDIM\,/PatchEdgeDataOpsIntegerNd\</g'                     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchEdgeDataOpsReal\<NDIM\,/PatchEdgeDataOpsRealNd\</g'                           `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

echo "  PatchFaceData"
perl -pi -e 's/PatchFaceDataBasicOps\<NDIM\,/PatchFaceDataBasicOpsNd\</g'                         `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchFaceDataMisFaceaneousOpsReal\<NDIM\,/PatchFaceDataMisFaceaneousOpsRealNd\</g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchFaceDataNormOpsComplex\<NDIM\,/PatchFaceDataNormOpsComplexNd\</g'             `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchFaceDataNormOpsReal\<NDIM\,/PatchFaceDataNormOpsRealNd\</g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchFaceDataOpsComplex\<NDIM\,/PatchFaceDataOpsComplexNd\</g'                     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchFaceDataOpsInteger\<NDIM\,/PatchFaceDataOpsIntegerNd\</g'                     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchFaceDataOpsReal\<NDIM\,/PatchFaceDataOpsRealNd\</g'                           `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

echo "  PatchNodeData"
perl -pi -e 's/PatchNodeDataBasicOps\<NDIM\,/PatchNodeDataBasicOpsNd\</g'                         `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchNodeDataMisNodeaneousOpsReal\<NDIM\,/PatchNodeDataMisNodeaneousOpsRealNd\</g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchNodeDataNormOpsComplex\<NDIM\,/PatchNodeDataNormOpsComplexNd\</g'             `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchNodeDataNormOpsReal\<NDIM\,/PatchNodeDataNormOpsRealNd\</g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchNodeDataOpsComplex\<NDIM\,/PatchNodeDataOpsComplexNd\</g'                     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchNodeDataOpsInteger\<NDIM\,/PatchNodeDataOpsIntegerNd\</g'                     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchNodeDataOpsReal\<NDIM\,/PatchNodeDataOpsRealNd\</g'                           `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

echo "  PatchSideData"
perl -pi -e 's/PatchSideDataBasicOps\<NDIM\,/PatchSideDataBasicOpsNd\</g'                         `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchSideDataMisSideaneousOpsReal\<NDIM\,/PatchSideDataMisSideaneousOpsRealNd\</g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchSideDataNormOpsComplex\<NDIM\,/PatchSideDataNormOpsComplexNd\</g'             `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchSideDataNormOpsReal\<NDIM\,/PatchSideDataNormOpsRealNd\</g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchSideDataOpsComplex\<NDIM\,/PatchSideDataOpsComplexNd\</g'                     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchSideDataOpsInteger\<NDIM\,/PatchSideDataOpsIntegerNd\</g'                     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/PatchSideDataOpsReal\<NDIM\,/PatchSideDataOpsRealNd\</g'                           `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

# mesh
echo "mesh"
perl -pi -e 's/BergerRigoutsos\<NDIM\>/BergerRigoutsosNd/g'                       `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/BoxGeneratorStrategy\<NDIM\>/BoxGeneratorStrategyNd/g'             `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/GriddingAlgorithm\<NDIM\>/GriddingAlgorithmNd/g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/LoadBalancer\<NDIM\>/LoadBalancerNd/g'                             `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/StandardTagAndInitialize\<NDIM\>/StandardTagAndInitializeNd/g'     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/StandardTagAndInitStrategy\<NDIM\>/StandardTagAndInitStrategyNd/g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/TagAndInitializeStrategy\<NDIM\>/TagAndInitializeStrategyNd/g'     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

# pdat
echo "pdat"
perl -pi -e 's/ArrayData\<NDIM\,/ArrayDataNd\</g'                              `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CellData\<NDIM\,/CellDataNd\</g'                                `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CellDataFactory\<NDIM\,/CellDataFactoryNd\</g'                  `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CellDoubleConstantRefine\<NDIM\>/CellDoubleConstantRefineNd/g'  `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CellGeometry\<NDIM\>/CellGeometryNd/g'                          `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CellIterator\<NDIM\>/CellIteratorNd/g'                          `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CellOverlap\<NDIM\>/CellOverlapNd/g'                            `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CellVariable\<NDIM\,/CellVariableNd\</g'                        `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/EdgeData\<NDIM\,/EdgeDataNd\</g'                                `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/EdgeDataFactory\<NDIM\,/EdgeDataFactoryNd\</g'                  `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/EdgeGeometry\<NDIM\>/EdgeGeometryNd/g'                          `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/EdgeIterator\<NDIM\>/EdgeIteratorNd/g'                          `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/EdgeOverlap\<NDIM\>/EdgeOverlapNd/g'                            `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/EdgeVariable\<NDIM\,/EdgeVariableNd\</g'                        `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/FaceData\<NDIM\,/FaceDataNd\</g'                                `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/FaceDataFactory\<NDIM\,/FaceDataFactoryNd\</g'                  `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/FaceGeometry\<NDIM\>/FaceGeometryNd/g'                          `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/FaceIterator\<NDIM\>/FaceIteratorNd/g'                          `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/FaceOverlap\<NDIM\>/FaceOverlapNd/g'                            `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/FaceVariable\<NDIM\,/FaceVariableNd\</g'                        `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/IndexData\<NDIM\,/IndexDataNd\</g'                              `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v LSetData.cpp)`
perl -pi -e 's/IndexDataFactory\<NDIM\,/IndexDataFactoryNd\</g'                `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v LSetData.cpp)`
perl -pi -e 's/IndexDataNode\<NDIM\,/IndexDataNodeNd\</g'                      `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v LSetData.cpp)`
perl -pi -e 's/IndexIterator\<NDIM\,/IndexIteratorNd\</g'                      `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v LSetData.cpp)`
perl -pi -e 's/IndexVariable\<NDIM\,/IndexVariableNd\</g'                      `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v LSetData.cpp)`
perl -pi -e 's/NodeData\<NDIM\,/NodeDataNd\</g'                                `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/NodeDataFactory\<NDIM\,/NodeDataFactoryNd\</g'                  `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/NodeGeometry\<NDIM\>/NodeGeometryNd/g'                          `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/NodeIterator\<NDIM\>/NodeIteratorNd/g'                          `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/NodeOverlap\<NDIM\>/NodeOverlapNd/g'                            `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/NodeVariable\<NDIM\,/NodeVariableNd\</g'                        `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuteredgeData\<NDIM\,/OuteredgeDataNd\</g'                      `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuteredgeDataFactory\<NDIM\,/OuteredgeDataFactoryNd\</g'        `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuteredgeGeometry\<NDIM\>/OuteredgeGeometryNd/g'                `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuteredgeIterator\<NDIM\>/OuteredgeIteratorNd/g'                `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuteredgeVariable\<NDIM\,/OuteredgeVariableNd\</g'              `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuterfaceData\<NDIM\,/OuterfaceDataNd\</g'                      `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuterfaceDataFactory\<NDIM\,/OuterfaceDataFactoryNd\</g'        `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuterfaceGeometry\<NDIM\>/OuterfaceGeometryNd/g'                `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuterfaceIterator\<NDIM\>/OuterfaceIteratorNd/g'                `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuterfaceVariable\<NDIM\,/OuterfaceVariableNd\</g'              `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuternodeData\<NDIM\,/OuternodeDataNd\</g'                      `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuternodeDataFactory\<NDIM\,/OuternodeDataFactoryNd\</g'        `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuternodeGeometry\<NDIM\>/OuternodeGeometryNd/g'                `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuternodeIterator\<NDIM\>/OuternodeIteratorNd/g'                `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OuternodeVariable\<NDIM\,/OuternodeVariableNd\</g'              `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OutersideData\<NDIM\,/OutersideDataNd\</g'                      `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OutersideDataFactory\<NDIM\,/OutersideDataFactoryNd\</g'        `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OutersideGeometry\<NDIM\>/OutersideGeometryNd/g'                `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OutersideIterator\<NDIM\>/OutersideIteratorNd/g'                `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/OutersideVariable\<NDIM\,/OutersideVariableNd\</g'              `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/SideData\<NDIM\,/SideDataNd\</g'                                `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/SideDataFactory\<NDIM\,/SideDataFactoryNd\</g'                  `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/SideGeometry\<NDIM\>/SideGeometryNd/g'                          `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/SideIterator\<NDIM\>/SideIteratorNd/g'                          `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/SideOverlap\<NDIM\>/SideOverlapNd/g'                            `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/SideVariable\<NDIM\,/SideVariableNd\</g'                        `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

# solv
echo "solv"
perl -pi -e 's/LocationIndexRobinBcCoefs\<NDIM\>/LocationIndexRobinBcCoefsNd/g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/RobinBcCoefStrategy\<NDIM\>/RobinBcCoefStrategyNd/g'             `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/SAMRAIVectorReal\<NDIM\,/SAMRAIVectorRealNd\</g'                 `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

# xfer
echo "xfer"
perl -pi -e 's/BoxGeometryFillPattern\<NDIM\>/BoxGeometryFillPatternNd/g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CoarsenAlgorithm\<NDIM\>/CoarsenAlgorithmNd/g'             `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CoarsenOperator\<NDIM\>/CoarsenOperatorNd/g'               `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CoarsenPatchStrategy\<NDIM\>/CoarsenPatchStrategyNd/g'     `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/CoarsenSchedule\<NDIM\>/CoarsenScheduleNd/g'               `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/Geometry\<NDIM\>/GeometryNd/g'                             `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/RefineAlgorithm\<NDIM\>/RefineAlgorithmNd/g'               `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/RefineClasses\<NDIM\>/RefineClassesNd/g'                   `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/RefineOperator\<NDIM\>/RefineOperatorNd/g'                 `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/RefinePatchStrategy\<NDIM\>/RefinePatchStrategyNd/g'       `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/RefineSchedule\<NDIM\>/RefineScheduleNd/g'                 `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`
perl -pi -e 's/VariableFillPattern\<NDIM\>/VariableFillPatternNd/g'       `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h)`

# this is a little lazy...
perl -pi -e 's/SAMRAI\:\:tbox\:\:ConstPointer\</IBTK\:\:SAMRAIConstPointer\</g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v contrib)`
perl -pi -e 's/tbox\:\:ConstPointer\</SAMRAIConstPointer\</g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v contrib)`
perl -pi -e 's/ConstPointer\</SAMRAIConstPointer\</g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v contrib)`
perl -pi -e 's/SAMRAISAMRAIConstPointer/SAMRAIConstPointer/g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v contrib)`

perl -pi -e 's/SAMRAI\:\:tbox\:\:Pointer\</IBTK\:\:SAMRAIPointer\</g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v contrib)`
perl -pi -e 's/tbox\:\:Pointer\</SAMRAIPointer\</g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v contrib)`
perl -pi -e 's/Pointer\</SAMRAIPointer\</g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v contrib)`
perl -pi -e 's/SAMRAISAMRAIPointer/SAMRAIPointer/g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v contrib)`
perl -pi -e 's/SAMRAISAMRAIConstPointer/SAMRAIConstPointer/g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v contrib)`
perl -pi -e 's/SAMRAIConstSAMRAIPointer/SAMRAIConstPointer/g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep -v contrib)`

perl -pi -e 's/IBTK::SAMRAIPointer/SAMRAIPointer/g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep ibtk)`
perl -pi -e 's/IBTK::SAMRAIConstPointer/SAMRAIConstPointer/g' `(find . -name \*.cpp -o -name \*.h | grep -v samrai_compatability.h | grep ibtk)`

# do not use "new" to initialize a SAMRAIPointer
perl -i~ -0pe 's/IBTK::SAMRAIPointer\<(.*)\>\s+([^,\s]+)\s+=\s+new\s+([^(]+)\(/IBTK::SAMRAIPointer<$1> $2 = IBTK::make_samrai_shared\<$3\>\(/gm' `(find . -name \*.h | grep -v samrai_compatability.h)`
perl -i~ -0pe 's/SAMRAIPointer\<(.*)\>\s+([^,\s]+)\s+=\s+new\s+([^(]+)\(/SAMRAIPointer<$1> $2 = make_samrai_shared\<$3\>\(/gm' `(find . -name \*.h | grep -v samrai_compatability.h)`

# use auto if the types on the LHS and the RHS are the same:
perl -i~ -0pe 's/IBTK::SAMRAIPointer\<(.*)\>\s+([^,\s]+)\s+=\s+new\s+\1\(/auto $2 = IBTK::make_samrai_shared\<$1\>\(/gm' `(find . -name \*.cpp | grep -v samrai_compatability.h)`
perl -i~ -0pe 's/SAMRAIPointer\<(.*)\>\s+([^,\s]+)\s+=\s+new\s+\1\(/auto $2 = make_samrai_shared\<$1\>\(/gm' `(find . -name \*.cpp | grep -v samrai_compatability.h)`

# but don't use auto if the types on the LHS and RHS are different; conversion of Pointer<T> is not robust:
perl -i~ -0pe 's/IBTK::SAMRAIPointer\<(.*)\>\s+([^,\s]+)\s+=\s+new\s+([^(]+)\(/IBTK::SAMRAIPointer\<$1\> $2 = IBTK::make_samrai_shared\<$3\>\(/gm' `(find . -name \*.cpp | grep -v samrai_compatability.h)`
perl -i~ -0pe 's/SAMRAIPointer\<(.*)\>\s+([^,\s]+)\s+=\s+new\s+([^(]+)\(/SAMRAIPointer\<$1\> $2 = make_samrai_shared\<$3\>\(/gm' `(find . -name \*.cpp | grep -v samrai_compatability.h)`
