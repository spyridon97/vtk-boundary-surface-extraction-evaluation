#include <vtkm/Version.h>
#include <vtkm/cont/CellSetPermutation.h>
#include <vtkm/cont/CellSetSingleType.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/filter/clean_grid/CleanGrid.h>
#include <vtkm/filter/geometry_refinement/Tetrahedralize.h>

#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkSMPThreadLocalObject.h>
#include <vtkSMPTools.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersionFull.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkmlib/ArrayConverters.h>
#include <vtkmlib/CellSetConverters.h>
#include <vtkmlib/DataArrayConverters.h>
#include <vtkmlib/DataSetConverters.h>
#include <vtkmlib/UnstructuredGridConverter.h>
#include <vtksys/SystemInformation.hxx>

#include "ExternalFacesHashCountFnv1a.h"
#include "ExternalFacesHashCountMinPointId.h"
#include "ExternalFacesHashFightFnv1a.h"
#include "ExternalFacesHashFightMinPointId.h"
#include "ExternalFacesHashSortFnv1a.h"
#include "ExternalFacesHashSortMinPointId.h"

#include "vtkDataSetSurfaceFilterSHash.h"
#include "vtkGeometryFilterPClassifier.h"
#include "vtkGeometryFilterPHash.h"
#include "vtkGeometryFilterSClassifier.h"

#include "Arguments.h"
#include "YamlWriter.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <random>
#include <sstream>
#include <vector>

auto ReadDataSet(const std::string& filename) -> vtkSmartPointer<vtkUnstructuredGrid>
{
  vtkNew<vtkXMLUnstructuredGridReader> reader;
  reader->SetFileName(filename.c_str());
  reader->Update();
  return reader->GetOutput();
}

auto RandomizeDataSet(const vtkSmartPointer<vtkUnstructuredGrid>& ug, YamlWriter& log)
  -> vtkSmartPointer<vtkUnstructuredGrid>
{
  // Create a list of indices to shuffle (representing point IDs)
  std::vector<vtkIdType> pointMap(ug->GetNumberOfPoints());
  std::iota(pointMap.begin(), pointMap.end(), 0); // Initialize to {0, 1, 2, ..., numPoints-1}

  // Use a random engine and shuffle the point IDs with the given seed
  const auto seed = static_cast<vtkm::UInt32>(std::time(nullptr));
  log.AddDictionaryEntry("randomize-seed", seed);

  // Shuffle the point IDs
  std::mt19937 randomEngine(seed);
  std::shuffle(pointMap.begin(), pointMap.end(), randomEngine);

  vtkPoints* inputPoints = ug->GetPoints();
  // Create a copy of the points
  vtkNew<vtkPoints> randomPoints;
  randomPoints->SetDataType(inputPoints->GetDataType());
  randomPoints->SetNumberOfPoints(inputPoints->GetNumberOfPoints());

  // update the points in the randomized order
  vtkSMPTools::For(0, inputPoints->GetNumberOfPoints(),
    [&](vtkIdType begin, vtkIdType end)
    {
      double point[3];
      for (vtkIdType index = begin; index < end; index++)
      {
        inputPoints->GetPoint(pointMap[index], point);
        randomPoints->SetPoint(index, point);
      }
    });

  vtkCellArray* inputCells = ug->GetCells();
  // randomize the cells
  auto offsets = vtk::TakeSmartPointer(
    vtkDataArray::CreateDataArray(inputCells->GetOffsetsArray()->GetDataType()));
  auto connectivity = vtk::TakeSmartPointer(
    vtkDataArray::CreateDataArray(inputCells->GetConnectivityArray()->GetDataType()));
  offsets->ShallowCopy(inputCells->GetOffsetsArray());
  connectivity->SetNumberOfTuples(inputCells->GetConnectivityArray()->GetNumberOfTuples());
  vtkNew<vtkCellArray> randomCells;
  randomCells->SetData(offsets, connectivity);

  // update the cells with the new point IDs
  vtkSMPThreadLocalObject<vtkIdList> tlPointIds;
  vtkSMPTools::For(0, inputCells->GetNumberOfCells(),
    [&](vtkIdType begin, vtkIdType end)
    {
      vtkIdList* pointIds = tlPointIds.Local();
      for (vtkIdType index = begin; index < end; index++)
      {
        inputCells->GetCellAtId(index, pointIds);
        for (vtkIdType i = 0; i < pointIds->GetNumberOfIds(); i++)
        {
          pointIds->SetId(i, pointMap[pointIds->GetId(i)]);
        }
        randomCells->ReplaceCellAtId(index, pointIds);
      }
    });

  // copy cell types
  vtkNew<vtkUnsignedCharArray> cellTypes;
  cellTypes->ShallowCopy(ug->GetCellTypesArray());

  auto randomUG = vtkSmartPointer<vtkUnstructuredGrid>::New();
  randomUG->SetPoints(randomPoints);
  randomUG->SetCells(cellTypes, randomCells);
  randomUG->GetCellData()->PassData(ug->GetCellData());

  return randomUG;
}

template <typename ExternalFacesAlgorithm>
auto RunVTKTrial(ExternalFacesAlgorithm* externalFaces, vtkUnstructuredGrid* inData,
  YamlWriter& log, bool firstRun = false) -> vtkm::Float64
{
  vtkm::cont::Timer timer;
  timer.Start();
  externalFaces->SetInputData(inData);
  externalFaces->Modified();
  externalFaces->Update();
  auto outData = externalFaces->GetOutput();
  timer.Stop();
  vtkm::Float64 elapsedTime = timer.GetElapsedTime();
  if (firstRun)
  {
    log.AddDictionaryEntry("num-output-points", outData->GetNumberOfPoints());
    log.AddDictionaryEntry("num-output-cells", outData->GetNumberOfCells());
  }
  return elapsedTime;
}

template <typename ExternalFacesAlgorithm>
auto DoVTKRun(const std::string& algorithmName, const std::string& hashName, unsigned int numTrials,
  vtkUnstructuredGrid* inData, YamlWriter& log) -> void
{
  vtkNew<ExternalFacesAlgorithm> externalFaces;
  log.StartListItem();
  log.AddDictionaryEntry("algorithm-name", algorithmName);
  log.AddDictionaryEntry("hash-name", hashName);
  log.AddDictionaryEntry("full-name", algorithmName + " " + hashName);

  log.AddDictionaryEntry(
    "first-run-time", RunVTKTrial(externalFaces.GetPointer(), inData, log, true));

  if (numTrials > 0)
  {
    log.StartBlock("trials");
    for (unsigned int trial = 0; trial < numTrials; trial++)
    {
      log.StartListItem();
      log.AddDictionaryEntry("trial-index", trial);
      log.AddDictionaryEntry("seconds-total", RunVTKTrial(externalFaces.GetPointer(), inData, log));
    }
    log.EndBlock();
  }
}

template <typename ExternalFacesWorklet>
auto RunVTKmTrial(ExternalFacesWorklet externalFaces, const vtkm::cont::DataSet& inData,
  YamlWriter& log, bool firstRun = false) -> vtkm::Float64
{
  const vtkm::cont::UnknownCellSet& unknownCellSet = inData.GetCellSet();
  auto inCellSet = unknownCellSet.ResetCellSetList<VTKM_DEFAULT_CELL_SET_LIST_UNSTRUCTURED>();

  vtkm::cont::CellSetExplicit<> outCellSet;

  std::stringstream dummyStream;
  YamlWriter dummyLog(dummyStream);

  vtkm::cont::Timer timer;
  timer.Start();
  externalFaces.Run(inCellSet, outCellSet, firstRun ? dummyLog : log);
  timer.Stop();
  vtkm::Float64 elapsedTime = timer.GetElapsedTime();
  vtkm::filter::clean_grid::CleanGrid cleanGrid;
  cleanGrid.SetMergePoints(false);
  cleanGrid.SetCompactPointFields(true);
  vtkm::cont::DataSet outDataSet;
  outDataSet.AddCoordinateSystem(inData.GetCoordinateSystem());
  outDataSet.SetCellSet(outCellSet);
  timer.Start();
  auto cleanResult = cleanGrid.Execute(outDataSet);
  timer.Stop();
  elapsedTime += timer.GetElapsedTime();
  if (firstRun)
  {
    log.AddDictionaryEntry(
      "num-output-points", cleanResult.GetCoordinateSystem().GetNumberOfPoints());
    log.AddDictionaryEntry("num-output-cells", cleanResult.GetNumberOfCells());
  }
  else
  {
    log.AddDictionaryEntry("seconds-clean-grid", timer.GetElapsedTime());
  }
  return elapsedTime;
}

template <typename ExternalFacesWorklet>
auto DoVTKmRun(const std::string& algorithmName, const std::string& hashName,
  unsigned int numTrials, const vtkm::cont::DataSet& inData, YamlWriter& log) -> void
{
  ExternalFacesWorklet externalFaces;
  log.StartListItem();
  log.AddDictionaryEntry("algorithm-name", algorithmName);
  log.AddDictionaryEntry("hash-name", hashName);
  log.AddDictionaryEntry("full-name", algorithmName + " " + hashName);

  log.AddDictionaryEntry("first-run-time", RunVTKmTrial(externalFaces, inData, log, true));

  if (numTrials > 0)
  {
    log.StartBlock("trials");
    for (unsigned int trial = 0; trial < numTrials; trial++)
    {
      log.StartListItem();
      log.AddDictionaryEntry("trial-index", trial);
      log.AddDictionaryEntry("seconds-total", RunVTKmTrial(externalFaces, inData, log));
    }
    log.EndBlock();
  }
}

auto ComputeFaceHashDistribution(
  const vtkSmartPointer<vtkUnstructuredGrid>& inData, YamlWriter& log) -> void
{
  const vtkIdType numPoints = inData->GetNumberOfPoints();
  std::vector<unsigned int> minPointCounter(inData->GetNumberOfPoints(), 0);
  std::vector<unsigned int> fnv1aCounter(inData->GetNumberOfPoints(), 0);

  vtkNew<vtkGenericCell> cell;
  for (vtkIdType cellId = 0; cellId < inData->GetNumberOfCells(); cellId++)
  {
    inData->GetCell(cellId, cell);
    for (int faceId = 0; faceId < cell->GetNumberOfFaces(); ++faceId)
    {
      vtkIdList* pointIds = cell->GetFace(faceId)->GetPointIds();
      std::sort(pointIds->GetPointer(0), pointIds->GetPointer(0) + pointIds->GetNumberOfIds());
      vtkIdType minPointId = pointIds->GetId(0);
      ++minPointCounter[minPointId];
      vtkm::Id3 canonicalFaceId(pointIds->GetId(0), pointIds->GetId(1), pointIds->GetId(2));
      vtkm::Id fnv1aHash = vtkm::Hash(canonicalFaceId) % numPoints;
      ++fnv1aCounter[fnv1aHash];
    }
  }

  std::map<unsigned int, unsigned int> minPointHashSizeMap;
  std::map<unsigned int, unsigned int> fnv1aHashSizeMap;
  for (vtkIdType pointId = 0; pointId < numPoints; pointId++)
  {
    auto minPointHashSizeMapIt = minPointHashSizeMap.find(minPointCounter[pointId]);
    if (minPointHashSizeMapIt == minPointHashSizeMap.end())
    {
      minPointHashSizeMap[minPointCounter[pointId]] = 1;
    }
    else
    {
      minPointHashSizeMapIt->second++;
    }
    if (fnv1aHashSizeMap.find(fnv1aCounter[pointId]) == fnv1aHashSizeMap.end())
    {
      fnv1aHashSizeMap[fnv1aCounter[pointId]] = 1;
    }
    else
    {
      fnv1aHashSizeMap[fnv1aCounter[pointId]]++;
    }
  }

  log.StartBlock("face-hash-distribution");
  log.StartBlock("MinPointId");
  for (const auto& entry : minPointHashSizeMap)
  {
    log.AddDictionaryEntry(std::to_string(entry.first), entry.second);
  }
  log.EndBlock();
  log.StartBlock("FNV1A");
  for (const auto& entry : fnv1aHashSizeMap)
  {
    log.AddDictionaryEntry(std::to_string(entry.first), entry.second);
  }
  log.EndBlock();
  log.EndBlock();
}

auto main(int argc, char** argv) -> int
{
  Arguments args;
  args.ParseArguments(argc, argv);

  vtksys::SystemInformation sysinfo;

  std::string deviceName =
    args.DeviceName != "TBB" ? vtkm::cont::make_DeviceAdapterId(args.DeviceName).GetName() : "TBB";

  YamlWriter log;
  log.StartListItem();

  log.AddDictionaryEntry("vtk-version", VTK_VERSION_FULL);
  log.AddDictionaryEntry("vtkm-version", VTKM_VERSION_FULL);
  log.AddDictionaryEntry("hostname", sysinfo.GetHostname());
  std::time_t currentTime = std::time(nullptr);
  char timeString[256];
  std::strftime(timeString, 256, "%Y-%m-%dT%H:%M:%S%z", std::localtime(&currentTime));
  log.AddDictionaryEntry("date", timeString);

  vtkSMPTools::Initialize(static_cast<int>(args.NumberOfThreads));
  // Construct the command line string for vtkm::cont::Initialize
  std::vector<std::string> strings = { argv[0], "--vtkm-device", deviceName };
  strings.emplace_back("--vtkm-num-threads");
  strings.push_back(std::to_string(args.NumberOfThreads));
  std::vector<char*> argvVector;
  for (const auto& str : strings)
  {
    argvVector.push_back(const_cast<char*>(str.c_str()));
  }
  argvVector.push_back(nullptr);
  int vtkm_argc = static_cast<int>(argvVector.size() - 1);
  char** vtkm_argv = argvVector.data();
  auto result = vtkm::cont::Initialize(vtkm_argc, vtkm_argv,
    vtkm::cont::InitializeOptions::RequireDevice | vtkm::cont::InitializeOptions::ErrorOnBadOption);
  log.AddDictionaryEntry("device", result.Device.GetName());
  log.AddDictionaryEntry("num-threads", args.NumberOfThreads);

  log.AddDictionaryEntry("input-file", args.InputFileName);

  vtkSmartPointer<vtkUnstructuredGrid> vtkInputData = ReadDataSet(args.InputFileName);

  if (args.Randomize)
  {
    vtkInputData = RandomizeDataSet(vtkInputData, log);
    log.AddDictionaryEntry("topology-connections", "randomized");
  }
  else
  {
    log.AddDictionaryEntry("topology-connections", "regular");
  }
  log.AddDictionaryEntry("num-input-points", vtkInputData->GetNumberOfPoints());
  log.AddDictionaryEntry("num-input-cells", vtkInputData->GetNumberOfCells());

  // Convert the VTK data to VTK-m data if needed
  // vtkm::cont::DataSet vtkmInputData;
  // if (args.PHashSort || args.PHashFight || args.PHashCount)
  // {
  //   vtkmInputData = tovtkm::Convert(vtkInputData, tovtkm::FieldsFlag::PointsAndCells);
  // }
  vtkm::cont::DataSet vtkmInputData =
    tovtkm::Convert(vtkInputData, tovtkm::FieldsFlag::PointsAndCells);
  // deallocate the VTK data if it is not needed
  // if (!(args.HashDistribution || args.SClassifier || args.SHash || args.PClassifier ||
  // args.PHash))
  // {
  //   vtkInputData = nullptr;
  // }

  const auto datasetMemoryUsed = sysinfo.GetProcMemoryUsed();
  log.AddDictionaryEntry("dataset-memory-used", datasetMemoryUsed);

  log.StartBlock("experiments");

  if (args.HashDistribution)
  {
    ComputeFaceHashDistribution(vtkInputData, log);
  }
  if (args.SClassifier)
  {
    DoVTKRun<vtkGeometryFilterSClassifier>(
      "S-Classifier", "None", args.NumberOfTrials, vtkInputData, log);
  }
  if (args.SHash)
  {
    DoVTKRun<vtkDataSetSurfaceFilterSHash>(
      "S-Hash", "MinPointID", args.NumberOfTrials, vtkInputData, log);
  }
  if (args.PClassifier)
  {
    DoVTKRun<vtkGeometryFilterPClassifier>(
      "P-Classifier", "None", args.NumberOfTrials, vtkInputData, log);
  }
  if (args.PHash)
  {
    DoVTKRun<vtkGeometryFilterPHash>(
      "P-Hash", "MinPointID", args.NumberOfTrials, vtkInputData, log);
  }

  if (args.PHashSort)
  {
    if (args.HashFunction == 0 || args.HashFunction == 1)
    {
      DoVTKmRun<vtkm::worklet::ExternalFacesHashSortFnv1a>(
        "P-Hash-Sort", "FNV1A", args.NumberOfTrials, vtkmInputData, log);
    }
    if (args.HashFunction == 0 || args.HashFunction == 2)
    {
      DoVTKmRun<vtkm::worklet::ExternalFacesHashSortMinPointId>(
        "P-Hash-Sort", "MinPointID", args.NumberOfTrials, vtkmInputData, log);
    }
  }
  if (args.PHashFight)
  {
    if (args.HashFunction == 0 || args.HashFunction == 1)
    {
      DoVTKmRun<vtkm::worklet::ExternalFacesHashFightFnv1a>(
        "P-Hash-Fight", "FNV1A", args.NumberOfTrials, vtkmInputData, log);
    }
    if (args.HashFunction == 0 || args.HashFunction == 2)
    {
      DoVTKmRun<vtkm::worklet::ExternalFacesHashFightMinPointId>(
        "P-Hash-Fight", "MinPointID", args.NumberOfTrials, vtkmInputData, log);
    }
  }
  if (args.PHashCount)
  {
    if (args.HashFunction == 0 || args.HashFunction == 1)
    {
      DoVTKmRun<vtkm::worklet::ExternalFacesHashCountFnv1a>(
        "P-Hash-Count", "FNV1A", args.NumberOfTrials, vtkmInputData, log);
    }
    if (args.HashFunction == 0 || args.HashFunction == 2)
    {
      DoVTKmRun<vtkm::worklet::ExternalFacesHashCountMinPointId>(
        "P-Hash-Count", "MinPointID", args.NumberOfTrials, vtkmInputData, log);
    }
  }
  log.EndBlock();
}