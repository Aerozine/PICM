# 🌊 PICM 🌊

A small 2D fluid simulation and method comparison:
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
detected automatically when available or fetched.

## Run One Case

```bash
./build-release/bin/PIC test/CH4/section-4-4-1/test-uniform.json
```

Each JSON file describes the grid, method, pressure solver, boundary
conditions, and output folder.
`test/CH*` mirrors the report chapters; `test/extra` keeps older runs that are not always rechecked.

## Run ALL the Report Cases

```bash
./run.sh
```
## Clean

```bash
make clean
```
