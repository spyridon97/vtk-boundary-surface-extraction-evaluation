//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_m_worklet_ExternalFacesHashCountFnv1a_h
#define vtk_m_worklet_ExternalFacesHashCountFnv1a_h

#include <vtkm/CellShape.h>
#include <vtkm/Hash.h>
#include <vtkm/Math.h>

#include <vtkm/exec/CellFace.h>

#include <vtkm/Swap.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include <vtkm/cont/CellSetExplicit.h>
#include <vtkm/cont/ConvertNumComponentsToOffsets.h>
#include <vtkm/cont/Timer.h>

#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/worklet/ScatterCounting.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletMapTopology.h>

#include "YamlWriter.h"

namespace vtkm
{
namespace worklet
{

struct ExternalFacesHashCountFnv1a
{
  // Worklet that returns the number of faces for each cell/shape
  class NumFacesPerCell : public vtkm::worklet::WorkletVisitCellsWithPoints
  {
  public:
    using ControlSignature = void(CellSetIn inCellSet, FieldOut facesInCell);
    using ExecutionSignature = void(CellShape, _2);
    using InputDomain = _1;

    template <typename CellShapeTag>
    VTKM_EXEC void operator()(CellShapeTag shape, vtkm::IdComponent& facesInCell) const
    {
      vtkm::exec::CellFaceNumberOfFaces(shape, facesInCell);
    }
  };

  // Worklet that identifies a cell face by a hash value. Not necessarily completely unique.
  class FaceHash : public vtkm::worklet::WorkletVisitCellsWithPoints
  {
  public:
    using ControlSignature = void(CellSetIn cellset, FieldOutCell cellFaceHashes);
    using ExecutionSignature = void(CellShape, PointIndices, _2);
    using InputDomain = _1;

    explicit FaceHash(const vtkm::Id& hashTableSize)
      : HashTableSize(hashTableSize)
    {
    }

    template <typename CellShapeTag, typename CellNodeVecType, typename CellFaceHashes>
    VTKM_EXEC void operator()(const CellShapeTag shape, const CellNodeVecType& cellNodeIds,
      CellFaceHashes& cellFaceHashes) const
    {
      const vtkm::IdComponent numFaces = cellFaceHashes.GetNumberOfComponents();
      for (vtkm::IdComponent faceIndex = 0; faceIndex < numFaces; ++faceIndex)
      {
        vtkm::Id3 faceId;
        vtkm::exec::CellFaceCanonicalId(faceIndex, shape, cellNodeIds, faceId);
        cellFaceHashes[faceIndex] = (vtkm::Hash(faceId) % this->HashTableSize);
      }
    }

  private:
    vtkm::Id HashTableSize;
  };

  // Worklet that identifies the number of faces per hash.
  class NumFacesPerHash : public vtkm::worklet::WorkletMapField
  {
  public:
    using ControlSignature = void(FieldIn cellFaceHashes, AtomicArrayInOut facesPerHash);
    using ExecutionSignature = void(_1, _2);
    using InputDomain = _1;

    template <typename CellFaceHashes, typename FacesPerHashArray>
    VTKM_EXEC void operator()(
      const CellFaceHashes& cellFaceHashes, FacesPerHashArray& facesPerHash) const
    {
      const vtkm::IdComponent numFaces = cellFaceHashes.GetNumberOfComponents();
      for (vtkm::IdComponent faceIndex = 0; faceIndex < numFaces; ++faceIndex)
      {
        // MemoryOrder::Relaxed is safe here, since we're not using the atomics for synchronization.
        facesPerHash.Add(cellFaceHashes[faceIndex], 1, vtkm::MemoryOrder::Relaxed);
      }
    }
  };

  // Worklet that writes out the cell and face ids of of each face per hash.
  class BuildFacesPerHash : public vtkm::worklet::WorkletMapField
  {
  public:
    using ControlSignature = void(FieldIn cellFaceHashes, AtomicArrayInOut facesPerHash,
      WholeArrayOut cellIdOfFacesPerHash, WholeArrayOut faceIdOfFacesPerHash);
    using ExecutionSignature = void(InputIndex, _1, _2, _3, _4);
    using InputDomain = _1;

    template <typename CellFaceHashes, typename FacesPerHashArray,
      typename CellIdOfFacePerHashArray, typename FaceIdOfFacePerHashArray>
    VTKM_EXEC void operator()(vtkm::Id inputIndex, const CellFaceHashes& cellFaceHashes,
      FacesPerHashArray& facesPerHash, CellIdOfFacePerHashArray& cellIdOfFacesPerHash,
      FaceIdOfFacePerHashArray& faceIdOfFacesPerHash) const
    {
      const vtkm::IdComponent numFaces = cellFaceHashes.GetNumberOfComponents();
      for (vtkm::IdComponent faceIndex = 0; faceIndex < numFaces; ++faceIndex)
      {
        const auto& faceHash = cellFaceHashes[faceIndex];
        // MemoryOrder::Relaxed is safe here, since we're not using the atomics for synchronization.
        const vtkm::IdComponent hashFaceIndex =
          facesPerHash.Add(faceHash, -1, vtkm::MemoryOrder::Relaxed) - 1;
        cellIdOfFacesPerHash.Get(faceHash)[hashFaceIndex] = inputIndex;
        faceIdOfFacesPerHash.Get(faceHash)[hashFaceIndex] = faceIndex;
      }
    }
  };

  // Worklet that identifies the number of external faces per Hash.
  // Because there can be collisions in the hash, this instance hash might
  // represent multiple faces, which have to be checked. The resulting
  // number is the total number of external faces. It also reorders the
  // faces so that the external faces are first, followed by the internal faces.
  class FaceCounts : public vtkm::worklet::WorkletMapField
  {
  public:
    using ControlSignature = void(FieldInOut cellIdOfFacesInHash, FieldInOut faceIdOfFacesInHash,
      WholeCellSetIn<> inputCells, FieldOut externalFacesInHash);
    using ExecutionSignature = void(_1, _2, _3, _4);
    using InputDomain = _1;

    template <typename CellIdOfFacesInHash, typename FaceIdOfFacesInHash, typename CellSetType>
    VTKM_EXEC void operator()(CellIdOfFacesInHash& cellIdOfFacesInHash,
      FaceIdOfFacesInHash& faceIdOfFacesInHash, const CellSetType& cellSet,
      vtkm::IdComponent& externalFacesInHash) const
    {
      using CellIdType = typename CellIdOfFacesInHash::ComponentType;
      using FaceIdType = typename FaceIdOfFacesInHash::ComponentType;

      const vtkm::IdComponent numFacesInHash = cellIdOfFacesInHash.GetNumberOfComponents();
      VTKM_ASSERT(faceIdOfFacesInHash.GetNumberOfComponents() == numFacesInHash);

      // Start by assuming all faces are duplicate, then remove two for each duplicate pair we find.
      externalFacesInHash = numFacesInHash;

      static constexpr vtkm::IdComponent FACE_CANONICAL_IDS_CACHE_SIZE = 100;
      if (numFacesInHash <= 1)
      {
        // Either one or zero faces. If there is one, it's external, In either case, do nothing.
      }
      else if (FACE_CANONICAL_IDS_CACHE_SIZE >= numFacesInHash) // fast path
      {
        vtkm::Vec<vtkm::Id3, FACE_CANONICAL_IDS_CACHE_SIZE> faceCanonicalIds;
        for (vtkm::IdComponent faceIndex = 0; faceIndex < numFacesInHash; ++faceIndex)
        {
          vtkm::exec::CellFaceCanonicalId(faceIdOfFacesInHash[faceIndex],
            cellSet.GetCellShape(cellIdOfFacesInHash[faceIndex]),
            cellSet.GetIndices(cellIdOfFacesInHash[faceIndex]), faceCanonicalIds[faceIndex]);
        }
        // iterate over the faces in the hash in reverse order (to minime the swaps being
        // performed), find duplicates faces, and put them at the end.
        vtkm::IdComponent newDuplicateIndex = numFacesInHash - 1;
        for (vtkm::IdComponent myIndex = newDuplicateIndex; myIndex >= 1;
             myIndex = vtkm::Min(myIndex - 1, newDuplicateIndex))
        {
          const auto& myFace = faceCanonicalIds[myIndex];
          for (vtkm::IdComponent otherIndex = myIndex - 1; otherIndex >= 0; --otherIndex)
          {
            const auto& otherFace = faceCanonicalIds[otherIndex];
            if (myFace == otherFace)
            {
              // Faces are the same. Must be internal. We don't have to worry about otherFace
              // matching anything else because a proper topology will have at most 2 cells sharing
              // a face, so there should be no more matches.
              externalFacesInHash -= 2;
              // move the duplicates/internal faces at the end to avoid revisiting them
              if (newDuplicateIndex != myIndex)
              {
                FaceCounts::SwapFace<CellIdType, FaceIdType>(cellIdOfFacesInHash[newDuplicateIndex],
                  faceIdOfFacesInHash[newDuplicateIndex], cellIdOfFacesInHash[myIndex],
                  faceIdOfFacesInHash[myIndex]);
                vtkm::Swap(faceCanonicalIds[newDuplicateIndex], faceCanonicalIds[myIndex]);
              }
              --newDuplicateIndex;
              if (newDuplicateIndex != otherIndex)
              {
                FaceCounts::SwapFace<CellIdType, FaceIdType>(cellIdOfFacesInHash[newDuplicateIndex],
                  faceIdOfFacesInHash[newDuplicateIndex], cellIdOfFacesInHash[otherIndex],
                  faceIdOfFacesInHash[otherIndex]);
                vtkm::Swap(faceCanonicalIds[newDuplicateIndex], faceCanonicalIds[otherIndex]);
              }
              --newDuplicateIndex;
              break;
            }
          }
        }
      }
      else // slow path
      {
        // iterate over the faces in the hash in reverse order (to minime the swaps being
        // performed), find duplicates faces, and put the them at the end.
        vtkm::IdComponent newDuplicateIndex = numFacesInHash - 1;
        for (vtkm::IdComponent myIndex = newDuplicateIndex; myIndex >= 1;
             myIndex = vtkm::Min(myIndex - 1, newDuplicateIndex))
        {
          vtkm::Id3 myFace;
          vtkm::exec::CellFaceCanonicalId(faceIdOfFacesInHash[myIndex],
            cellSet.GetCellShape(cellIdOfFacesInHash[myIndex]),
            cellSet.GetIndices(cellIdOfFacesInHash[myIndex]), myFace);

          for (vtkm::IdComponent otherIndex = myIndex - 1; otherIndex >= 0; --otherIndex)
          {
            vtkm::Id3 otherFace;
            vtkm::exec::CellFaceCanonicalId(faceIdOfFacesInHash[otherIndex],
              cellSet.GetCellShape(cellIdOfFacesInHash[otherIndex]),
              cellSet.GetIndices(cellIdOfFacesInHash[otherIndex]), otherFace);
            if (myFace == otherFace)
            {
              // Faces are the same. Must be internal. We don't have to worry about otherFace
              // matching anything else because a proper topology will have at most 2 cells sharing
              // a face, so there should be no more matches.
              externalFacesInHash -= 2;
              // move the duplicates/internal faces at the end to avoid revisiting them
              if (newDuplicateIndex != myIndex)
              {
                FaceCounts::SwapFace<CellIdType, FaceIdType>(cellIdOfFacesInHash[newDuplicateIndex],
                  faceIdOfFacesInHash[newDuplicateIndex], cellIdOfFacesInHash[myIndex],
                  faceIdOfFacesInHash[myIndex]);
              }
              --newDuplicateIndex;
              if (newDuplicateIndex != otherIndex)
              {
                FaceCounts::SwapFace<CellIdType, FaceIdType>(cellIdOfFacesInHash[newDuplicateIndex],
                  faceIdOfFacesInHash[newDuplicateIndex], cellIdOfFacesInHash[otherIndex],
                  faceIdOfFacesInHash[otherIndex]);
              }
              --newDuplicateIndex;
              break;
            }
          }
        }
      }
    }

    template <typename T1Type, typename T2Type, typename T1Reference, typename T2Reference>
    VTKM_EXEC inline static void SwapFace(
      T1Reference&& cell1, T2Reference&& face1, T1Reference&& cell2, T2Reference&& face2)
    {
      const T1Type tmpCell = cell1;
      cell1 = cell2;
      cell2 = tmpCell;
      const T2Type tmpFace = face1;
      face1 = face2;
      face2 = tmpFace;
    }
  };

public:
  // Worklet that returns the number of points for each outputted face.
  // Have to manage the case where multiple faces have the same hash.
  class NumPointsPerFace : public vtkm::worklet::WorkletMapField
  {
  public:
    using ControlSignature = void(FieldIn cellIdOfFacesInHash, FieldIn faceIdOfFacesInHash,
      WholeCellSetIn<> inputCells, FieldOut pointsInExternalFace);
    using ExecutionSignature = void(_1, _2, _3, VisitIndex, _4);
    using InputDomain = _1;

    using ScatterType = vtkm::worklet::ScatterCounting;

    template <typename CountArrayType>
    VTKM_CONT static ScatterType MakeScatter(const CountArrayType& countArray)
    {
      VTKM_IS_ARRAY_HANDLE(CountArrayType);
      return ScatterType(countArray);
    }

    template <typename CellIdOfFacesInHash, typename FaceIdOfFacesInHash, typename CellSetType>
    VTKM_EXEC void operator()(const CellIdOfFacesInHash& cellIdOfFacesInHash,
      const FaceIdOfFacesInHash& faceIdOfFacesInHash, const CellSetType& cellSet,
      vtkm::IdComponent visitIndex, vtkm::IdComponent& pointsInExternalFace) const
    {
      // external faces are first, so we can use the visit index directly
      vtkm::exec::CellFaceNumberOfPoints(faceIdOfFacesInHash[visitIndex],
        cellSet.GetCellShape(cellIdOfFacesInHash[visitIndex]), pointsInExternalFace);
    }
  };

  // Worklet that returns the shape and connectivity for each external face
  class BuildConnectivity : public vtkm::worklet::WorkletMapField
  {
  public:
    using ControlSignature = void(FieldIn cellIdOfFacesInHash, FieldIn faceIdOfFacesInHash,
      WholeCellSetIn<> inputCells, FieldOut shapesOut, FieldOut connectivityOut,
      FieldOut cellIdMapOut);
    using ExecutionSignature = void(_1, _2, _3, VisitIndex, _4, _5, _6);
    using InputDomain = _1;

    using ScatterType = vtkm::worklet::ScatterCounting;

    template <typename CellIdOfFacesInHash, typename FaceIdOfFacesInHash, typename CellSetType,
      typename ConnectivityType>
    VTKM_EXEC void operator()(const CellIdOfFacesInHash& cellIdOfFacesInHash,
      const FaceIdOfFacesInHash& faceIdOfFacesInHash, const CellSetType& cellSet,
      vtkm::IdComponent visitIndex, vtkm::UInt8& shapeOut, ConnectivityType& connectivityOut,
      vtkm::Id& cellIdMapOut) const
    {
      // external faces are first, so we can use the visit index directly
      const typename CellIdOfFacesInHash::ComponentType myCellId = cellIdOfFacesInHash[visitIndex];
      const typename FaceIdOfFacesInHash::ComponentType myFaceId = faceIdOfFacesInHash[visitIndex];

      const typename CellSetType::CellShapeTag shapeIn = cellSet.GetCellShape(myCellId);
      vtkm::exec::CellFaceShape(myFaceId, shapeIn, shapeOut);
      cellIdMapOut = myCellId;

      vtkm::IdComponent numFacePoints;
      vtkm::exec::CellFaceNumberOfPoints(myFaceId, shapeIn, numFacePoints);
      VTKM_ASSERT(numFacePoints == connectivityOut.GetNumberOfComponents());

      const typename CellSetType::IndicesType inCellIndices = cellSet.GetIndices(myCellId);
      for (vtkm::IdComponent facePointIndex = 0; facePointIndex < numFacePoints; ++facePointIndex)
      {
        vtkm::IdComponent localFaceIndex;
        const vtkm::ErrorCode status =
          vtkm::exec::CellFaceLocalIndex(facePointIndex, myFaceId, shapeIn, localFaceIndex);
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
  ExternalFacesHashCountFnv1a() {}

  void ReleaseCellMapArrays() { this->CellIdMap.ReleaseResources(); }

  ///////////////////////////////////////////////////
  /// \brief ExternalFacesHashCountFnv1a: Extract Faces on outside of geometry
  template <typename InCellSetType, typename ShapeStorage, typename ConnectivityStorage,
    typename OffsetsStorage>
  VTKM_CONT void Run(const InCellSetType& inCellSet,
    vtkm::cont::CellSetExplicit<ShapeStorage, ConnectivityStorage, OffsetsStorage>& outCellSet,
    YamlWriter& log)
  {
    using PointCountArrayType = vtkm::cont::ArrayHandle<vtkm::IdComponent>;
    using ShapeArrayType = vtkm::cont::ArrayHandle<vtkm::UInt8, ShapeStorage>;
    using OffsetsArrayType = vtkm::cont::ArrayHandle<vtkm::Id, OffsetsStorage>;
    using ConnectivityArrayType = vtkm::cont::ArrayHandle<vtkm::Id, ConnectivityStorage>;

    // create an array to store the number of faces for each cell
    vtkm::cont::ArrayHandle<vtkm::IdComponent> facesPerCell;

    // compute the number of faces for each cell
    vtkm::worklet::DispatcherMapTopology<NumFacesPerCell> numFacesDispatcher;
    vtkm::cont::Timer timer;
    timer.Start();
    numFacesDispatcher.Invoke(inCellSet, facesPerCell);
    timer.Stop();
    log.AddDictionaryEntry("seconds-num-faces-per-cell", timer.GetElapsedTime());

    // Compute the offsets of the prefix sum of the number of faces per cell
    vtkm::Id totalNumberOfFaces;
    vtkm::cont::ArrayHandle<vtkm::Id> facesPerCellOffsets;
    timer.Start();
    vtkm::cont::ConvertNumComponentsToOffsets(
      facesPerCell, facesPerCellOffsets, totalNumberOfFaces);
    timer.Stop();
    log.AddDictionaryEntry("seconds-face-per-cell-count", timer.GetElapsedTime());
    // Release the resources of facesPerCell that is not needed anymore
    facesPerCell.ReleaseResources();

    if (totalNumberOfFaces == 0)
    {
      // Data has no faces. Output is empty.
      outCellSet.PrepareToAddCells(0, 0);
      outCellSet.CompleteAddingCells(inCellSet.GetNumberOfPoints());
      return;
    }
    // allocate the array to store the hash values of the faces
    vtkm::cont::ArrayHandle<vtkm::HashType> faceHashes;
    faceHashes.Allocate(totalNumberOfFaces);

    // create a group vec array to access the faces of each cell conveniently
    auto faceHashesGroupVec =
      vtkm::cont::make_ArrayHandleGroupVecVariable(faceHashes, facesPerCellOffsets);

    // compute the hash values of the faces
    const vtkm::Id numberOfHashes = inCellSet.GetNumberOfPoints();
    vtkm::worklet::DispatcherMapTopology<FaceHash> faceHashDispatcher((FaceHash(numberOfHashes)));
    timer.Start();
    faceHashDispatcher.Invoke(inCellSet, faceHashesGroupVec);
    timer.Stop();
    log.AddDictionaryEntry("seconds-face-hash", timer.GetElapsedTime());

    // We create an ArrayHandle and pass it to the Worklet as AtomicArrayInOut.
    vtkm::cont::ArrayHandle<vtkm::IdComponent> facesPerHash;
    facesPerHash.AllocateAndFill(numberOfHashes, 0);

    // count the number of faces for each hash
    vtkm::worklet::DispatcherMapField<NumFacesPerHash> numFacesPerHashDispatcher;
    timer.Start();
    numFacesPerHashDispatcher.Invoke(faceHashesGroupVec, facesPerHash);
    timer.Stop();
    log.AddDictionaryEntry("seconds-num-faces-per-hash", timer.GetElapsedTime());

    // compute the offsets of the prefix sum of the number of faces per hash
    vtkm::cont::ArrayHandle<vtkm::Id> facesPerHashOffsets;
    timer.Start();
    vtkm::cont::ConvertNumComponentsToOffsets(facesPerHash, facesPerHashOffsets);
    timer.Stop();
    log.AddDictionaryEntry("seconds-face-per-hash-count", timer.GetElapsedTime());

    // we create the arrays to store the cell and face ids of each face per hash
    vtkm::cont::ArrayHandle<vtkm::Id> cellIdOfFacesPerHash;
    cellIdOfFacesPerHash.Allocate(totalNumberOfFaces);
    vtkm::cont::ArrayHandle<vtkm::Int8> faceIdOfFacesPerHash;
    faceIdOfFacesPerHash.Allocate(totalNumberOfFaces);

    // create a group vec array to access/write the cell ids of each face per hash
    auto cellIdOfFacesPerHashGroupVec =
      vtkm::cont::make_ArrayHandleGroupVecVariable(cellIdOfFacesPerHash, facesPerHashOffsets);
    // create a group vec array to access/write the face ids of each face per hash
    auto faceIdOfFacesPerHashGroupVec =
      vtkm::cont::make_ArrayHandleGroupVecVariable(faceIdOfFacesPerHash, facesPerHashOffsets);

    // Build the cell and face ids of all faces per hash
    vtkm::worklet::DispatcherMapField<BuildFacesPerHash> buildFacesPerHashDispatcher;
    timer.Start();
    buildFacesPerHashDispatcher.Invoke(
      faceHashesGroupVec, facesPerHash, cellIdOfFacesPerHashGroupVec, faceIdOfFacesPerHashGroupVec);
    timer.Stop();
    log.AddDictionaryEntry("seconds-build-faces-per-hash", timer.GetElapsedTime());
    // Release the resources of the arrays that are not needed anymore
    facesPerCellOffsets.ReleaseResources();
    faceHashes.ReleaseResources();
    facesPerHash.ReleaseResources();

    // create an array to count the number of external faces per hash
    vtkm::cont::ArrayHandle<vtkm::IdComponent> externalFacesPerHash;
    externalFacesPerHash.Allocate(numberOfHashes);

    // compute the number of external faces per hash
    vtkm::worklet::DispatcherMapField<FaceCounts> faceCountsDispatcher;
    timer.Start();
    faceCountsDispatcher.Invoke(
      cellIdOfFacesPerHashGroupVec, faceIdOfFacesPerHashGroupVec, inCellSet, externalFacesPerHash);
    timer.Stop();
    log.AddDictionaryEntry("seconds-face-counts", timer.GetElapsedTime());

    // create a scatter array to access the hashes with external faces
    timer.Start();
    auto scatterCullInternalFaces = NumPointsPerFace::MakeScatter(externalFacesPerHash);
    timer.Stop();
    log.AddDictionaryEntry("seconds-scatter-cull-internal-faces", timer.GetElapsedTime());
    const vtkm::Id numberOfExternalFaces = scatterCullInternalFaces.GetOutputRange(numberOfHashes);
    // Release the resources of externalFacesPerHash that is not needed anymore
    externalFacesPerHash.ReleaseResources();

    // create an array to store the number of points of the external faces
    PointCountArrayType pointsPerExternalFace;
    pointsPerExternalFace.Allocate(numberOfExternalFaces);

    // compute the number of points of the external faces
    vtkm::worklet::DispatcherMapField<NumPointsPerFace> pointsPerFaceDispatcher(
      scatterCullInternalFaces);
    timer.Start();
    pointsPerFaceDispatcher.Invoke(
      cellIdOfFacesPerHashGroupVec, faceIdOfFacesPerHashGroupVec, inCellSet, pointsPerExternalFace);
    timer.Stop();
    log.AddDictionaryEntry("seconds-points-per-face", timer.GetElapsedTime());

    // create an array to store the shape of the external faces
    ShapeArrayType externalFacesShapes;
    externalFacesShapes.Allocate(numberOfExternalFaces);

    // compute the offsets of the prefix sum of the number of points per eternal face
    OffsetsArrayType pointsPerExternalFaceOffsets;
    vtkm::Id connectivitySize;
    timer.Start();
    vtkm::cont::ConvertNumComponentsToOffsets(
      pointsPerExternalFace, pointsPerExternalFaceOffsets, connectivitySize);
    timer.Stop();
    log.AddDictionaryEntry("seconds-face-point-count", timer.GetElapsedTime());

    // create an array to connectivity of the external faces
    ConnectivityArrayType externalFacesConnectivity;
    externalFacesConnectivity.Allocate(connectivitySize);

    // create a group vec array to access the connectivity of each external face
    auto externalFacesConnectivityGroupVec = vtkm::cont::make_ArrayHandleGroupVecVariable(
      externalFacesConnectivity, pointsPerExternalFaceOffsets);

    // create an array to store the cell id of the external faces
    vtkm::cont::ArrayHandle<vtkm::Id> faceToCellIdMap;
    faceToCellIdMap.Allocate(numberOfExternalFaces);

    // build the connectivity of the external faces
    vtkm::worklet::DispatcherMapField<BuildConnectivity> buildConnectivityDispatcher(
      scatterCullInternalFaces);
    timer.Start();
    buildConnectivityDispatcher.Invoke(cellIdOfFacesPerHashGroupVec, faceIdOfFacesPerHashGroupVec,
      inCellSet, externalFacesShapes, externalFacesConnectivityGroupVec, faceToCellIdMap);
    timer.Stop();
    log.AddDictionaryEntry("seconds-build-connectivity", timer.GetElapsedTime());

    outCellSet.Fill(inCellSet.GetNumberOfPoints(), externalFacesShapes, externalFacesConnectivity,
      pointsPerExternalFaceOffsets);
    this->CellIdMap = faceToCellIdMap;
  }

  vtkm::cont::ArrayHandle<vtkm::Id> GetCellIdMap() const { return this->CellIdMap; }

private:
  vtkm::cont::ArrayHandle<vtkm::Id> CellIdMap;

}; // struct ExternalFacesHashCountFnv1a
};
} // namespace vtkm::worklet

#endif // vtk_m_worklet_ExternalFacesHashCountFnv1a_h
