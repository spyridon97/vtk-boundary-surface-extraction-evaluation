# VTK Boundary Surface Extraction Evaluation

Repository to evaluate the performance of boundary surface extraction filters in VTK

## Build Instructions

There are 5 directories in this repository which test different implementations of the boundary surface extraction
filter.

* seqClassifier
* parClassifier
* seqHashing
* parHashing
* parGroupedHashing

### seqClassifier

```bash
cd seqClassifier
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DVTK_SMP_IMPLEMENTATION_TYPE=TBB ..
cmake --build . --target seqClassifier -j
``` 

### parClassifier

```bash
cd parClassifier
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DVTK_SMP_IMPLEMENTATION_TYPE=TBB ..
cmake --build . --target parClassifier -j
```

### seqHashing

```bash
cd seqHashing
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DVTK_SMP_IMPLEMENTATION_TYPE=TBB ..
cmake --build . --target seqHashing -j
```

### parHashing

```bash
cd parHashing
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DVTK_SMP_IMPLEMENTATION_TYPE=TBB ..
cmake --build . --target parHashing -j
```

### parGroupedHashing

```bash
cd parGroupedHashing
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DVTK_SMP_IMPLEMENTATION_TYPE=TBB ..
cmake --build . --target parGroupedHashing -j
```

## Datasets

https://drive.google.com/drive/folders/12TfEi9djQM3TTk8LkpERdQjxLQXONXCp?usp=sharing

## Execution Instructions

### seqClassifier

```bash
/usr/bin/time -vv ./seqClassifier/build/seqClassifier inputFile numberOfIterations  
```

### parClassifier

```bash
/usr/bin/time -vv ./parClassifier/build/parClassifier inputFile numberOfIterations numberOfThreads  
```

### seqHashing

```bash
/usr/bin/time -vv ./seqHashing/build/seqHashing inputFile numberOfIterations  
```

### parHashing

```bash
/usr/bin/time -vv ./parHashing/build/parHashing inputFile numberOfIterations numberOfThreads  
```

### parGroupedHashing

```bash
/usr/bin/time -vv ./parGroupedHashing/build/parGroupedHashing inputFile numberOfIterations numberOfThreads  
```
