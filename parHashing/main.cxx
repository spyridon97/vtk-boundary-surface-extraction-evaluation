/*=========================================================================

  Program:   Visualization Toolkit
  Module:    main.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <vtkNew.h>
#include <vtkSMPTools.h>
#include <vtkGeometryFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtksys/SystemInformation.hxx>

#include <chrono>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
  if (argc < 4)
  {
    std::cerr << "Usage: " << argv[0] << " filename numberOfIterations numberOfThreads"
              << std::endl;
    return EXIT_FAILURE;
  }
  const std::string filename = argv[1];
  int numberOfIterations = std::stoi(argv[2]);
  int numberOfThreads = std::stoi(argv[3]);

  if (numberOfThreads <= 1)
  {
    vtkSMPTools::SetBackend("SEQUENTIAL");
    vtkSMPTools::Initialize(1);
  }
  else
  {
    vtkSMPTools::SetBackend("TBB");
    vtkSMPTools::Initialize(numberOfThreads);
  }

  vtkNew<vtkXMLUnstructuredGridReader> reader;
  reader->SetFileName(filename.c_str());
  reader->Update();

  auto output = reader->GetOutput();

  std::cout << "Number of input cells: " << output->GetNumberOfCells() << std::endl;
  std::cout << "Number of input points: " << output->GetNumberOfPoints() << std::endl;

  vtksys::SystemInformation sysinfo;
  const auto memoryUsedBeforeFilter = sysinfo.GetProcMemoryUsed();
  std::cout << "Memory used by dataset in KB: " << memoryUsedBeforeFilter << std::endl;
  size_t totalTime = 0;
  for (int i = 0; i < numberOfIterations; ++i)
  {
    auto start = std::chrono::high_resolution_clock::now();
    auto filter = vtkSmartPointer<vtkGeometryFilter>::New();
    filter->SetInputConnection(reader->GetOutputPort());
    filter->Update();
    auto end = std::chrono::high_resolution_clock::now();
    totalTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Number of output cells: " << filter->GetOutput()->GetNumberOfCells() << std::endl;
    filter = nullptr;
  }
  std::cout << "Time: " << totalTime / numberOfIterations << std::endl;

  return EXIT_SUCCESS;
}
