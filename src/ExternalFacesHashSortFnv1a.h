//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_m_worklet_ExternalFacesHashSortFnv1a_h
#define vtk_m_worklet_ExternalFacesHashSortFnv1a_h

#include <vtkm/CellShape.h>
#include <vtkm/Hash.h>
#include <vtkm/Math.h>

#include <vtkm/exec/CellFace.h>

#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayCopyDevice.h>
#include <vtkm/cont/ArrayGetValues.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleConcatenate.h>
#include <vtkm/cont/ArrayHandleConstant.h>
#include <vtkm/cont/ArrayHandleGroupVec.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/ArrayHandlePermutation.h>
#include <vtkm/cont/ArrayHandleTransform.h>
#include <vtkm/cont/ArrayHandleView.h>
#include <vtkm/cont/CellSetExplicit.h>
#include <vtkm/cont/ConvertNumComponentsToOffsets.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/Field.h>
#include <vtkm/cont/Timer.h>

#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/worklet/DispatcherReduceByKey.h>
#include <vtkm/worklet/Keys.h>
#include <vtkm/worklet/ScatterCounting.h>
#include <vtkm/worklet/WorkletMapTopology.h>
#include <vtkm/worklet/WorkletReduceByKey.h>

#include "YamlWriter.h"

namespace vtkm
{
namespace worklet
{

struct ExternalFacesHashSortFnv1a
{
  // Worklet that returns the number of faces for each cell/shape
  class NumFacesPerCell : public vtkm::worklet::WorkletVisitCellsWithPoints
  {
  public:
    using ControlSignature = void(CellSetIn inCellSet, FieldOut numFacesInCell);
    using ExecutionSignature = void(CellShape, _2);
    using InputDomain = _1;

    template <typename CellShapeTag>
    VTKM_EXEC void operator()(CellShapeTag shape, vtkm::IdComponent& numFaces) const
    {
      vtkm::exec::CellFaceNumberOfFaces(shape, numFaces);
    }
  };

  // Worklet that identifies a cell face by a hash value. Not necessarily completely unique.
  class FaceHash : public vtkm::worklet::WorkletVisitCellsWithPoints
  {
  public:
    using ControlSignature = void(CellSetIn cellset,
                                  FieldOut faceHashes,
                                  FieldOut originCells,
                                  FieldOut originFaces);
    using ExecutionSignature = void(_2, _3, _4, CellShape, PointIndices, InputIndex, VisitIndex);
    using InputDomain = _1;

    using ScatterType = vtkm::worklet::ScatterCounting;

    template <typename CellShapeTag, typename CellNodeVecType>
    VTKM_EXEC void operator()(vtkm::HashType& faceHash,
                              vtkm::Id& cellIndex,
                              vtkm::IdComponent& faceIndex,
                              CellShapeTag shape,
                              const CellNodeVecType& cellNodeIds,
                              vtkm::Id inputIndex,
                              vtkm::IdComponent visitIndex) const
    {
      vtkm::Id3 faceId;
      vtkm::exec::CellFaceCanonicalId(visitIndex, shape, cellNodeIds, faceId);
      faceHash = vtkm::Hash(faceId);

      cellIndex = inputIndex;
      faceIndex = visitIndex;
    }
  };

  // Worklet that identifies the number of cells written out per face.
  // Because there can be collisions in the face ids, this instance might
  // represent multiple faces, which have to be checked. The resulting
  // number is the total number of external faces.
  class FaceCounts : public vtkm::worklet::WorkletReduceByKey
  {
  public:
    using ControlSignature = void(KeysIn keys,
                                  WholeCellSetIn<> inputCells,
                                  ValuesIn originCells,
                                  ValuesIn originFaces,
                                  ReducedValuesOut numOutputCells);
    using ExecutionSignature = _5(_2, _3, _4);
    using InputDomain = _1;

    template <typename CellSetType, typename OriginCellsType, typename OriginFacesType>
    VTKM_EXEC vtkm::IdComponent operator()(const CellSetType& cellSet,
                                           const OriginCellsType& originCells,
                                           const OriginFacesType& originFaces) const
    {
      vtkm::IdComponent numCellsOnHash = originCells.GetNumberOfComponents();
      VTKM_ASSERT(originFaces.GetNumberOfComponents() == numCellsOnHash);

      // Start by assuming all faces are unique, then remove one for each
      // face we find a duplicate for.
      vtkm::IdComponent numExternalFaces = numCellsOnHash;

      for (vtkm::IdComponent myIndex = 0;
           myIndex < numCellsOnHash - 1; // Don't need to check last face
           myIndex++)
      {
        vtkm::Id3 myFace;
        vtkm::exec::CellFaceCanonicalId(originFaces[myIndex],
                                        cellSet.GetCellShape(originCells[myIndex]),
                                        cellSet.GetIndices(originCells[myIndex]),
                                        myFace);
        for (vtkm::IdComponent otherIndex = myIndex + 1; otherIndex < numCellsOnHash; otherIndex++)
        {
          vtkm::Id3 otherFace;
          vtkm::exec::CellFaceCanonicalId(originFaces[otherIndex],
                                          cellSet.GetCellShape(originCells[otherIndex]),
                                          cellSet.GetIndices(originCells[otherIndex]),
                                          otherFace);
          if (myFace == otherFace)
          {
            // Faces are the same. Must be internal. Remove 2, one for each face. We don't have to
            // worry about otherFace matching anything else because a proper topology will have at
            // most 2 cells sharing a face, so there should be no more matches.
            numExternalFaces -= 2;
            break;
          }
        }
      }

      return numExternalFaces;
    }
  };

private:
  // Resolves duplicate hashes by finding a specified unique face for a given hash.
  // Given a cell set (from a WholeCellSetIn) and the cell/face id pairs for each face
  // associated with a given hash, returns the index of the cell/face provided of the
  // visitIndex-th unique face. Basically, this method searches through all the cell/face
  // pairs looking for unique sets and returns the one associated with visitIndex.
  template <typename CellSetType, typename OriginCellsType, typename OriginFacesType>
  VTKM_EXEC static vtkm::IdComponent FindUniqueFace(const CellSetType& cellSet,
                                                    const OriginCellsType& originCells,
                                                    const OriginFacesType& originFaces,
                                                    vtkm::IdComponent visitIndex)
  {
    vtkm::IdComponent numCellsOnHash = originCells.GetNumberOfComponents();
    VTKM_ASSERT(originFaces.GetNumberOfComponents() == numCellsOnHash);

    // Find the visitIndex-th unique face.
    vtkm::IdComponent numFound = 0;
    vtkm::IdComponent myIndex = 0;
    while (true)
    {
      VTKM_ASSERT(myIndex < numCellsOnHash);
      vtkm::Id3 myFace;
      vtkm::exec::CellFaceCanonicalId(originFaces[myIndex],
                                      cellSet.GetCellShape(originCells[myIndex]),
                                      cellSet.GetIndices(originCells[myIndex]),
                                      myFace);
      bool foundPair = false;
      for (vtkm::IdComponent otherIndex = 0; otherIndex < numCellsOnHash; otherIndex++)
      {
        if (otherIndex == myIndex)
        {
          continue;
        }
        vtkm::Id3 otherFace;
        vtkm::exec::CellFaceCanonicalId(originFaces[otherIndex],
                                        cellSet.GetCellShape(originCells[otherIndex]),
                                        cellSet.GetIndices(originCells[otherIndex]),
                                        otherFace);
        if (myFace == otherFace)
        {
          // Faces are the same. Must be internal.
          foundPair = true;
          break;
        }
      }

      if (!foundPair)
      {
        if (numFound == visitIndex)
        {
          break;
        }
        else
        {
          numFound++;
        }
      }

      myIndex++;
    }

    return myIndex;
  }

public:
  // Worklet that returns the number of points for each outputted face.
  // Have to manage the case where multiple faces have the same hash.
  class NumPointsPerFace : public vtkm::worklet::WorkletReduceByKey
  {
  public:
    using ControlSignature = void(KeysIn keys,
                                  WholeCellSetIn<> inputCells,
                                  ValuesIn originCells,
                                  ValuesIn originFaces,
                                  ReducedValuesOut numPointsInFace);
    using ExecutionSignature = void(_2, _3, _4, VisitIndex, _5);
    using InputDomain = _1;

    using ScatterType = vtkm::worklet::ScatterCounting;

    template <typename CountArrayType>
    VTKM_CONT static ScatterType MakeScatter(const CountArrayType& countArray)
    {
      VTKM_IS_ARRAY_HANDLE(CountArrayType);
      return ScatterType(countArray);
    }

    template <typename CellSetType, typename OriginCellsType, typename OriginFacesType>
    VTKM_EXEC void operator()(const CellSetType& cellSet,
                              const OriginCellsType& originCells,
                              const OriginFacesType& originFaces,
                              vtkm::IdComponent visitIndex,
                              vtkm::IdComponent& numFacePoints) const
    {
      vtkm::IdComponent myIndex =
        ExternalFacesHashSortFnv1a::FindUniqueFace(cellSet, originCells, originFaces, visitIndex);

      vtkm::exec::CellFaceNumberOfPoints(
        originFaces[myIndex], cellSet.GetCellShape(originCells[myIndex]), numFacePoints);
    }
  };

  // Worklet that returns the shape and connectivity for each external face
  class BuildConnectivity : public vtkm::worklet::WorkletReduceByKey
  {
  public:
    using ControlSignature = void(KeysIn keys,
                                  WholeCellSetIn<> inputCells,
                                  ValuesIn originCells,
                                  ValuesIn originFaces,
                                  ReducedValuesOut shapesOut,
                                  ReducedValuesOut connectivityOut,
                                  ReducedValuesOut cellIdMapOut);
    using ExecutionSignature = void(_2, _3, _4, VisitIndex, _5, _6, _7);
    using InputDomain = _1;

    using ScatterType = vtkm::worklet::ScatterCounting;

    template <typename CellSetType,
              typename OriginCellsType,
              typename OriginFacesType,
              typename ConnectivityType>
    VTKM_EXEC void operator()(const CellSetType& cellSet,
                              const OriginCellsType& originCells,
                              const OriginFacesType& originFaces,
                              vtkm::IdComponent visitIndex,
                              vtkm::UInt8& shapeOut,
                              ConnectivityType& connectivityOut,
                              vtkm::Id& cellIdMapOut) const
    {
      const vtkm::IdComponent myIndex =
        ExternalFacesHashSortFnv1a::FindUniqueFace(cellSet, originCells, originFaces, visitIndex);
      const vtkm::IdComponent myFace = originFaces[myIndex];

      typename CellSetType::CellShapeTag shapeIn = cellSet.GetCellShape(originCells[myIndex]);
      vtkm::exec::CellFaceShape(myFace, shapeIn, shapeOut);
      cellIdMapOut = originCells[myIndex];

      vtkm::IdComponent numFacePoints;
      vtkm::exec::CellFaceNumberOfPoints(myFace, shapeIn, numFacePoints);

      VTKM_ASSERT(numFacePoints == connectivityOut.GetNumberOfComponents());

      typename CellSetType::IndicesType inCellIndices = cellSet.GetIndices(originCells[myIndex]);

      for (vtkm::IdComponent facePointIndex = 0; facePointIndex < numFacePoints; facePointIndex++)
      {
        vtkm::IdComponent localFaceIndex;
        vtkm::ErrorCode status =
          vtkm::exec::CellFaceLocalIndex(facePointIndex, myFace, shapeIn, localFaceIndex);
        if (status == vtkm::ErrorCode::Success)
        {
          connectivityOut[facePointIndex] = inCellIndices[localFaceIndex];
        }
        else
        {
          // An error condition, but do we want to crash the operation?
          connectivityOut[facePointIndex] = 0;
        }
      }
    }
  };

public:
  VTKM_CONT
  ExternalFacesHashSortFnv1a() {}

  void ReleaseCellMapArrays() { this->CellIdMap.ReleaseResources(); }

  ///////////////////////////////////////////////////
  /// \brief ExternalFaces: Extract Faces on outside of geometry
  template <typename InCellSetType,
            typename ShapeStorage,
            typename ConnectivityStorage,
            typename OffsetsStorage>
  VTKM_CONT void Run(
    const InCellSetType& inCellSet,
    vtkm::cont::CellSetExplicit<ShapeStorage, ConnectivityStorage, OffsetsStorage>& outCellSet,
    YamlWriter& log)
  {
    using PointCountArrayType = vtkm::cont::ArrayHandle<vtkm::IdComponent>;
    using ShapeArrayType = vtkm::cont::ArrayHandle<vtkm::UInt8, ShapeStorage>;
    using OffsetsArrayType = vtkm::cont::ArrayHandle<vtkm::Id, OffsetsStorage>;
    using ConnectivityArrayType = vtkm::cont::ArrayHandle<vtkm::Id, ConnectivityStorage>;

    // Create a worklet to map the number of faces to each cell
    vtkm::cont::ArrayHandle<vtkm::IdComponent> facesPerCell;
    vtkm::worklet::DispatcherMapTopology<NumFacesPerCell> numFacesDispatcher;

    vtkm::cont::Timer timer;
    timer.Start();
    numFacesDispatcher.Invoke(inCellSet, facesPerCell);
    timer.Stop();
    log.AddDictionaryEntry("seconds-num-faces-per-cell", timer.GetElapsedTime());

    timer.Start();
    vtkm::worklet::ScatterCounting scatterCellToFace(facesPerCell);
    timer.Stop();
    log.AddDictionaryEntry("seconds-face-input-count", timer.GetElapsedTime());
    facesPerCell.ReleaseResources();

    if (scatterCellToFace.GetOutputRange(inCellSet.GetNumberOfCells()) == 0)
    {
      // Data has no faces. Output is empty.
      outCellSet.PrepareToAddCells(0, 0);
      outCellSet.CompleteAddingCells(inCellSet.GetNumberOfPoints());
      return;
    }

    vtkm::cont::ArrayHandle<vtkm::HashType> faceHashes;
    vtkm::cont::ArrayHandle<vtkm::Id> originCells;
    vtkm::cont::ArrayHandle<vtkm::IdComponent> originFaces;
    vtkm::worklet::DispatcherMapTopology<FaceHash> faceHashDispatcher(scatterCellToFace);

    timer.Start();
    faceHashDispatcher.Invoke(inCellSet, faceHashes, originCells, originFaces);
    timer.Stop();
    log.AddDictionaryEntry("seconds-face-hash", timer.GetElapsedTime());

    timer.Start();
    vtkm::worklet::Keys<vtkm::HashType> faceKeys(faceHashes);
    timer.Stop();
    log.AddDictionaryEntry("seconds-keys-build-arrays", timer.GetElapsedTime());

    vtkm::cont::ArrayHandle<vtkm::IdComponent> faceOutputCount;
    vtkm::worklet::DispatcherReduceByKey<FaceCounts> faceCountDispatcher;

    timer.Start();
    faceCountDispatcher.Invoke(faceKeys, inCellSet, originCells, originFaces, faceOutputCount);
    timer.Stop();
    log.AddDictionaryEntry("seconds-face-count", timer.GetElapsedTime());

    timer.Start();
    auto scatterCullInternalFaces = NumPointsPerFace::MakeScatter(faceOutputCount);
    timer.Stop();
    log.AddDictionaryEntry("seconds-face-output-count", timer.GetElapsedTime());

    PointCountArrayType facePointCount;
    vtkm::worklet::DispatcherReduceByKey<NumPointsPerFace> pointsPerFaceDispatcher(
      scatterCullInternalFaces);

    timer.Start();
    pointsPerFaceDispatcher.Invoke(faceKeys, inCellSet, originCells, originFaces, facePointCount);
    timer.Stop();
    log.AddDictionaryEntry("seconds-points-per-face", timer.GetElapsedTime());

    ShapeArrayType faceShapes;

    OffsetsArrayType faceOffsets;
    vtkm::Id connectivitySize;
    timer.Start();
    vtkm::cont::ConvertNumComponentsToOffsets(facePointCount, faceOffsets, connectivitySize);
    timer.Stop();
    log.AddDictionaryEntry("seconds-face-point-count", timer.GetElapsedTime());

    ConnectivityArrayType faceConnectivity;
    // Must pre allocate because worklet invocation will not have enough
    // information to.
    faceConnectivity.Allocate(connectivitySize);

    vtkm::worklet::DispatcherReduceByKey<BuildConnectivity> buildConnectivityDispatcher(
      scatterCullInternalFaces);

    vtkm::cont::ArrayHandle<vtkm::Id> faceToCellIdMap;

    timer.Start();
    buildConnectivityDispatcher.Invoke(
      faceKeys,
      inCellSet,
      originCells,
      originFaces,
      faceShapes,
      vtkm::cont::make_ArrayHandleGroupVecVariable(faceConnectivity, faceOffsets),
      faceToCellIdMap);
    timer.Stop();
    log.AddDictionaryEntry("seconds-build-connectivity", timer.GetElapsedTime());

    outCellSet.Fill(inCellSet.GetNumberOfPoints(), faceShapes, faceConnectivity, faceOffsets);
    this->CellIdMap = faceToCellIdMap;
  }

  vtkm::cont::ArrayHandle<vtkm::Id> GetCellIdMap() const { return this->CellIdMap; }

private:
  vtkm::cont::ArrayHandle<vtkm::Id> CellIdMap;

}; // struct ExternalFaces
}
} // namespace vtkm::worklet

#endif // vtk_m_worklet_ExternalFacesHashSortFnv1a_h
