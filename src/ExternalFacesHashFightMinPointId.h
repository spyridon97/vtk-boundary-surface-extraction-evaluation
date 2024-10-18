//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_m_worklet_ExternalFacesHashFightMinPointId_h
#define vtk_m_worklet_ExternalFacesHashFightMinPointId_h

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

#include "vtkm/worklet/WorkletMapField.h"
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/worklet/DispatcherReduceByKey.h>
#include <vtkm/worklet/Keys.h>
#include <vtkm/worklet/ScatterCounting.h>
#include <vtkm/worklet/WorkletMapTopology.h>
#include <vtkm/worklet/WorkletReduceByKey.h>

#include "CellFaceMinPointId.h"
#include "YamlWriter.h"

namespace vtkm
{
namespace worklet
{

struct ExternalFacesHashFightMinPointId
{
  //Unary predicate operator
  //Returns True if the argument is equal to the constructor
  //integer argument; False otherwise.
  struct IsIntValue
  {
  private:
    int Value;

  public:
    VTKM_EXEC_CONT
    IsIntValue(const int& v)
      : Value(v)
    {
    }

    template <typename T>
    VTKM_EXEC_CONT bool operator()(const T& x) const
    {
      return x == T(Value);
    }
  };

  //Worklet that returns the number of faces for each cell/shape
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

  //Worklet that identifies a cell face by a hash value. Not necessarily completely unique.
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
      vtkm::Id minFacePointId;
      vtkm::exec::CellFaceMinPointId(visitIndex, shape, cellNodeIds, minFacePointId);
      faceHash = static_cast<vtkm::HashType>(minFacePointId);

      cellIndex = inputIndex;
      faceIndex = visitIndex;
    }
  };

  //Worklet that writes the face index at the location of the hash table.
  //Multiple entries are likely to write to the hash table, so they fight
  //and (hopefully) one wins.
  class HashFight : public vtkm::worklet::WorkletMapField
  {
  public:
    using ControlSignature = void(FieldIn Hashes, FieldIn FaceIds, WholeArrayOut HashTable);
    typedef void ExecutionSignature(_1, _2, _3);

    VTKM_CONT
    HashFight(vtkm::Id hashTableSize)
      : HashTableSize(hashTableSize)
    {
    }

    template <typename HashTablePortalType>
    VTKM_EXEC void operator()(vtkm::Id hash,
                              vtkm::Id faceId,
                              HashTablePortalType& hashTablePortal) const
    {
      hashTablePortal.Set(hash % this->HashTableSize, faceId);
    }

  private:
    vtkm::Id HashTableSize;
  };

  //Worklet that detects whether a face is internal.  If the
  //face is internal, then a value should not be assigned to the
  //face in the output array handle of face vertices; only external
  //faces should have a vector not equal to <-1,-1,-1>
  class CheckForMatches : public vtkm::worklet::WorkletMapField
  {
  public:
    typedef void ControlSignature(FieldIn activeHashes,
                                  FieldIn activeFaceIndices,
                                  WholeCellSetIn<> cellSet,
                                  WholeArrayIn originCells,
                                  WholeArrayIn originFaces,
                                  WholeArrayIn hashTable,
                                  FieldInOut isInactive,
                                  WholeArrayInOut isExternalFace);
    typedef void ExecutionSignature(_1, _2, _3, _4, _5, _6, _7, _8);

    VTKM_CONT
    CheckForMatches(vtkm::Id hashTableSize)
      : HashTableSize(hashTableSize)
    {
    }

    template <typename CellSetType,
              typename OriginCellsPortal,
              typename OriginFacesPortal,
              typename HashTablePortal,
              typename IsExternalFacePortal>
    VTKM_EXEC void operator()(vtkm::UInt32 hash,
                              vtkm::Id faceIndex,
                              const CellSetType& cellSet,
                              const OriginCellsPortal& originCellsPortal,
                              const OriginFacesPortal& originFacesPortal,
                              const HashTablePortal& hashTablePortal,
                              vtkm::UInt8& isInactive,
                              IsExternalFacePortal& isExternalFacePortal) const
    {
      vtkm::Id hashWinnerFace = hashTablePortal.Get(hash % this->HashTableSize);

      if (hashWinnerFace == faceIndex)
      {
        // Case 1: I won the hash fight by writing my index. I'm done so mark
        // myself as inactive.
        isInactive = vtkm::UInt8(1);
      }
      else
      {
        // Get a cononical representation of my face.
        vtkm::Id myOriginCell = originCellsPortal.Get(faceIndex);
        vtkm::IdComponent myOriginFace = originFacesPortal.Get(faceIndex);
        vtkm::Id3 myFace;
        vtkm::exec::CellFaceCanonicalId(myOriginFace,
                                        cellSet.GetCellShape(myOriginCell),
                                        cellSet.GetIndices(myOriginCell),
                                        myFace);

        // Get a cononical representation of the face in the hash table.
        vtkm::Id otherOriginCell = originCellsPortal.Get(hashWinnerFace);
        vtkm::IdComponent otherOriginFace = originFacesPortal.Get(hashWinnerFace);
        vtkm::Id3 otherFace;
        vtkm::exec::CellFaceCanonicalId(otherOriginFace,
                                        cellSet.GetCellShape(otherOriginCell),
                                        cellSet.GetIndices(otherOriginCell),
                                        otherFace);

        // See if these are the same face
        if (myFace == otherFace)
        {
          // Case 2: The faces are the same. This must be an internal face.
          // Mark both myself and the other face as internal.
          isInactive = vtkm::UInt8(1);
          isExternalFacePortal.Set(faceIndex, vtkm::UInt8(0));
          isExternalFacePortal.Set(hashWinnerFace, vtkm::UInt8(0));
        }
        else
        {
          // Case 3: I didn't win and my face didn't match. I didn't learn
          // anything so do nothing.
        }
      }
    }

  private:
    vtkm::Id HashTableSize;
  };

  // Worklet that counts the number of points that are in each (active) face.
  class NumPointsPerFace : public vtkm::worklet::WorkletMapField
  {
  public:
    typedef void ControlSignature(FieldIn faceIndices,
                                  WholeCellSetIn<> cellSet,
                                  WholeArrayIn originCells,
                                  WholeArrayIn originFaces,
                                  FieldOut numPointsInFace);
    typedef _5 ExecutionSignature(_1, _2, _3, _4);

    using ScatterType = vtkm::worklet::ScatterCounting;

    template <typename CellSetType, typename OriginCellsPortalType, typename OriginFacesPortalType>
    VTKM_EXEC vtkm::IdComponent operator()(vtkm::Id faceIndex,
                                           const CellSetType& cellSet,
                                           const OriginCellsPortalType& originCellsPortal,
                                           const OriginFacesPortalType& originFacesPortal) const
    {
      vtkm::Id originCell = originCellsPortal.Get(faceIndex);
      vtkm::IdComponent originFace = originFacesPortal.Get(faceIndex);

      vtkm::IdComponent numFacePoints;
      vtkm::exec::CellFaceNumberOfPoints(
        originFace, cellSet.GetCellShape(originCell), numFacePoints);
      return numFacePoints;
    }
  };

  // Worklet that writes out the shape and indices for each (active) face.
  class BuildConnectivity : public vtkm::worklet::WorkletMapField
  {
  public:
    typedef void ControlSignature(FieldIn faceIndices,
                                  WholeCellSetIn<> cellSet,
                                  WholeArrayIn originCells,
                                  WholeArrayIn originFaces,
                                  FieldOut shapesOut,
                                  FieldOut connectivityOut,
                                  FieldOut cellIdMapOut);
    typedef void ExecutionSignature(_1, _2, _3, _4, _5, _6, _7);

    using ScatterType = vtkm::worklet::ScatterCounting;

    template <typename CellSetType,
              typename OriginCellsPortalType,
              typename OriginFacesPortalType,
              typename ConnectivityType>
    VTKM_EXEC void operator()(const vtkm::Id faceIndex,
                              const CellSetType& cellSet,
                              const OriginCellsPortalType& originCellsPortal,
                              const OriginFacesPortalType& originFacesPortal,
                              vtkm::UInt8& shapeOut,
                              ConnectivityType& connectivityOut,
                              vtkm::Id& cellIdMapOut) const
    {
      vtkm::Id originCell = originCellsPortal.Get(faceIndex);
      vtkm::IdComponent originFace = originFacesPortal.Get(faceIndex);

      vtkm::exec::CellFaceShape(originFace, cellSet.GetCellShape(originCell), shapeOut);
      cellIdMapOut = originCell;

      vtkm::IdComponent numFacePoints;
      vtkm::exec::CellFaceNumberOfPoints(
        originFace, cellSet.GetCellShape(originCell), numFacePoints);
      VTKM_ASSERT(numFacePoints == connectivityOut.GetNumberOfComponents());

      typename CellSetType::IndicesType inCellIndices = cellSet.GetIndices(originCell);

      for (vtkm::IdComponent facePointIndex = 0; facePointIndex < numFacePoints; facePointIndex++)
      {
        vtkm::IdComponent localFaceIndex;
        vtkm::ErrorCode status = vtkm::exec::CellFaceLocalIndex(
          facePointIndex, originFace, cellSet.GetCellShape(originCell), localFaceIndex);
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
  ExternalFacesHashFightMinPointId() {}

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

    //Create a worklet to map the number of faces to each cell
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


    vtkm::Id totalNumFaces = faceHashes.GetNumberOfValues();

    timer.Start();

    //Set constant factor for hash table size: factor*totalFaces
    const vtkm::Id hashTableFactor = 2;

    vtkm::cont::ArrayHandle<vtkm::UInt8> isExternalFace;
    vtkm::cont::Algorithm::Copy(vtkm::cont::ArrayHandleConstant<vtkm::UInt8>(1, totalNumFaces),
                                isExternalFace);

    vtkm::cont::ArrayHandle<vtkm::Id> activeFaceIndices;
    vtkm::cont::Algorithm::Copy(vtkm::cont::ArrayHandleIndex(totalNumFaces), activeFaceIndices);

    vtkm::Id numActiveFaces = totalNumFaces;

    while (numActiveFaces > 0)
    {
      // Create a packe arrays of active face hashes
      auto activeHashes = vtkm::cont::make_ArrayHandlePermutation(activeFaceIndices, faceHashes);

      // Get ready the isInactive array.
      vtkm::cont::ArrayHandle<vtkm::UInt8> isInactive;
      vtkm::cont::Algorithm::Copy(vtkm::cont::ArrayHandleConstant<vtkm::UInt8>(0, numActiveFaces),
                                  isInactive);

      vtkm::Id hashTableSize = numActiveFaces * hashTableFactor;

      vtkm::cont::ArrayHandle<vtkm::Id> hashTable;
      hashTable.Allocate(hashTableSize);

      // Have all active hashes try to write their index to the hash table
      vtkm::worklet::DispatcherMapField<HashFight> fightDispatcher((HashFight(hashTableSize)));
      fightDispatcher.Invoke(activeHashes, activeFaceIndices, hashTable);


      // Have all active faces check to see if they matched and update
      // isInactive/isExternalFace.
      vtkm::worklet::DispatcherMapField<CheckForMatches> matchDispatcher(
        (CheckForMatches(hashTableSize)));
      matchDispatcher.Invoke(activeHashes,
                             activeFaceIndices,
                             inCellSet,
                             originCells,
                             originFaces,
                             hashTable,
                             isInactive,
                             isExternalFace);

      // Compact the activeFaceIndices by the isInactive flag.
      vtkm::cont::ArrayHandle<vtkm::Id> compactedActiveFaceIndices;
      vtkm::cont::Algorithm::CopyIf(
        activeFaceIndices, isInactive, compactedActiveFaceIndices, IsIntValue(0));
      activeFaceIndices = compactedActiveFaceIndices;

      // Update the number of active faces
      numActiveFaces = activeFaceIndices.GetNumberOfValues();
    }
    timer.Stop();
    log.AddDictionaryEntry("seconds-hash-fight-iterations", timer.GetElapsedTime());

    vtkm::worklet::ScatterCounting scatterCullInternalFaces(isExternalFace);

    PointCountArrayType facePointCount;
    vtkm::worklet::DispatcherMapField<NumPointsPerFace> pointsPerFaceDispatcher(
      scatterCullInternalFaces);

    timer.Start();
    pointsPerFaceDispatcher.Invoke(vtkm::cont::ArrayHandleIndex(totalNumFaces),
                                   inCellSet,
                                   originCells,
                                   originFaces,
                                   facePointCount);
    timer.Stop();
    log.AddDictionaryEntry("seconds-face-output-count", timer.GetElapsedTime());

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

    vtkm::worklet::DispatcherMapField<BuildConnectivity> buildConnectivityDispatcher(
      scatterCullInternalFaces);

    vtkm::cont::ArrayHandle<vtkm::Id> faceToCellIdMap;

    // Create a view that doesn't have the last offset:
    auto faceOffsetsTrim =
      vtkm::cont::make_ArrayHandleView(faceOffsets, 0, faceOffsets.GetNumberOfValues() - 1);

    timer.Start();
    buildConnectivityDispatcher.Invoke(
      vtkm::cont::ArrayHandleIndex(totalNumFaces),
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

}; //struct ExternalFacesHashFightMinPointId
}
} //namespace vtkm::worklet

#endif //vtk_m_worklet_ExternalFacesHashFightMinPointId_h
