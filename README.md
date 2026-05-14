# PICM

A small C++ codebase for 2D fluid simulation and method comparison:
Semi-Lagrangian, PIC, FLIP, and APIC.

Outputs are written as VTK files and can be opened directly in ParaView.

## Clone

```bash
git clone https://github.com/Aerozine/PICM.git
cd PICM
```

## Build

```bash
cmake -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-release
```

You need CMake, Ninja, and a C++17 compiler. Eigen, OpenMP, Zlib, and CUDA are
detected automatically when available.

## Run One Case

```bash
./build-release/bin/PIC test/CH4/section-4-4-1/test-uniform.json
```

Each JSON file describes the grid, method, pressure solver, boundary
conditions, and output folder.

## Run Report Cases

```bash
./run.sh
```

The script builds in release mode, then runs every `test/CH*/**/*.json` case,
following the same idea as the `heatsink` targets in the Makefile.

## Clean

```bash
make clean
```
