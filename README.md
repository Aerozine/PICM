# 🌊 PICM 🌊
[![CMake CI](https://github.com/Aerozine/PICM/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/Aerozine/PICM/actions/workflows/cmake-multi-platform.yml)

PICM is a 2D fluid simulation playground for comparing several grid and
particle methods:

- Semi-Lagrangian (`sl`)
- Vanilla PIC (`pic`)
- FLIP (`flip`)
- APIC (`apic`)

The codebase also includes several pressure solvers, optional OpenMP
parallelism, optional CUDA acceleration for selected solvers, and VTK output
for post-processing in ParaView.

## Current status

Implemented simulation methods:

- `sl`
- `pic`
- `flip`
- `apic`

Implemented pressure solvers:

- `jacobi`
- `gauss_seidel`
- `red_black_gauss_seidel`
- `cg`
- `miccg0`

Notes:

- `mixed_flip_pic` is declared in the code, but it is not implemented yet.
- CUDA acceleration is currently available for `red_black_gauss_seidel` and
  `cg`.

## Dependencies

Required:

- CMake 3.18 or newer
- A C++17 compiler

Optional:

- OpenMP for CPU parallelism
- Eigen3 for CPU `cg` and `miccg0`
- Zlib for compressed VTK output
- CUDA Toolkit for GPU pressure solvers

Dependency handling:

- `nlohmann/json` is fetched automatically by CMake if it is not available
  locally.
- Eigen3, OpenMP, Zlib, and CUDA are detected if present.

## Build

### CPU build

Release build:

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release -j
```

Debug build:

```bash
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug -j
```

Important debug behavior:

- If `CMAKE_BUILD_TYPE=Debug` and you do not override them explicitly,
  `USE_PARALLEL` and `USE_GPU` default to `OFF`.

### WSL / CUDA build

If `nvcc` is not on `PATH`, pass the compiler explicitly:

```bash
cmake -S . -B build-wsl \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_GPU=ON \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda-13.1/bin/nvcc \
  -DCUDAToolkit_ROOT=/usr/local/cuda-13.1

cmake --build build-wsl -j
```

You can also point CMake to another CUDA installation, for example an NVIDIA
HPC SDK `nvcc`, by changing `CMAKE_CUDA_COMPILER` and `CUDAToolkit_ROOT`.

## CMake options

Main configuration switches:

- `USE_FLOAT_PRECISION=ON|OFF`
- `USE_FAST_ALG=ON|OFF`
- `USE_PARALLEL=ON|OFF`
- `USE_GPU=ON|OFF`
- `GPU_ARCH=<number or sm_XX>`

Examples:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_PARALLEL=ON \
  -DUSE_GPU=OFF \
  -DUSE_FLOAT_PRECISION=ON
```

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_GPU=ON \
  -DGPU_ARCH=89
```

## Run

The executable expects the JSON config path as its only argument:

```bash
./build-release/bin/PIC test/PIC/downbreak.json
```

The legacy wrapper form is also accepted:

```bash
./build-release/bin/PIC -c test/PIC/downbreak.json
```

On Windows:

```powershell
.\build-release\bin\PIC.exe test\PIC\downbreak.json
```

The program prints the parsed configuration, runs the selected solver, and
writes outputs into the folder specified by the JSON file.

## Example cases

The repository already contains ready-to-run scenes under `test/`:

- `test/SL/`
- `test/PIC/`
- `test/FLIP/`
- `test/APIC/`

Some useful examples:

- `test/SL/small-von-karman.json`
- `test/SL/large-von-karman.json`
- `test/PIC/downbreak.json`
- `test/FLIP/freeFall.json`
- `test/APIC/downbreak.json`

## Configuration file

Simulation parameters are described in a JSON file. Common fields include:

- `dx`, `dy`, `dt`
- `nx`, `ny`, `nt`
- `density`
- `sampling_rate`
- `folder`
- `freeSurface`
- `gravity`
- `max_cfl`
- `ppcx`, `ppcy` for particle methods
- `coefPic` for FLIP blending
- `kernelOrder`
- output toggles such as `write_u`, `write_v`, `write_p`,
  `write_vorticity`, `write_particles`

### Method selection

The top-level `method` field selects the simulation family:

- `"sl"`
- `"pic"`
- `"flip"`
- `"apic"`

The nested `solver.type` field selects the pressure solver:

- `"jacobi"`
- `"gauss_seidel"`
- `"red_black_gauss_seidel"`
- `"cg"`
- `"miccg0"`

### Scene description

The scene can contain:

- `fluid`
- `air`
- `solid`
- `velocityu`
- `velocityv`
- `pressure`
- `smoke`

Supported object shapes:

- `rectangle`
- `cylinder`
- `u_tube` / `utube` / `u-tube` / `tube_u` / `manometer` / `manometre`

Coordinates can be integers or simple expressions such as `"nx/2"` or
`"ny-10"`.

`u_tube` draws a rounded U-shaped channel on the grid. `bottom_y` is the
lowest inner point of the bend; `left_x` and `right_x` place the two vertical
legs. Use the same geometry under `solid` to create the tube walls and under
`fluid` to initialise the liquid:

```json
"solid": {
  "u_tube": {
    "left_x": "55",
    "right_x": "155",
    "bottom_y": "14",
    "top_y": "165",
    "tube_width": "30",
    "wall": "6"
  }
},
"fluid": {
  "u_tube": {
    "left_x": "55",
    "right_x": "155",
    "bottom_y": "14",
    "top_y": "165",
    "tube_width": "30",
    "left_level": "150",
    "right_level": "112"
  }
}
```

### Minimal example

```json
{
  "dx": 0.05,
  "dy": 0.05,
  "dt": 0.01,
  "nx": 100,
  "ny": 120,
  "nt": 200,
  "density": 1000,
  "sampling_rate": 5,
  "folder": "results/PIC/example",
  "method": "pic",
  "freeSurface": true,
  "gravity": 9.81,
  "ppcx": 4,
  "ppcy": 4,
  "max_cfl": 0.75,
  "write_particles": true,
  "solver": {
    "type": "red_black_gauss_seidel",
    "max_iterations": 10000,
    "tolerance": 1e-3
  },
  "fluid": {
    "rectangle": {
      "x1": "1",
      "y1": "1",
      "x2": "nx/3",
      "y2": "ny/2"
    }
  },
  "solid": {
    "rectangle": [
      { "x1": "0", "y1": "0", "x2": "nx+1", "y2": "0" },
      { "x1": "0", "y1": "ny+1", "x2": "nx+1", "y2": "ny+1" },
      { "x1": "0", "y1": "0", "x2": "0", "y2": "ny+1" },
      { "x1": "nx+1", "y1": "0", "x2": "nx+1", "y2": "ny+1" }
    ]
  }
}
```

## Output

Depending on the enabled flags, the solver can write:

- `u`
- `v`
- `p`
- `div`
- `normVelocity`
- `vorticity`
- `smoke`
- `particles`
- `label`

Output format:

- grid fields are written as `.vti`
- particles are written as `.vtp`
- time collections are indexed through `.pvd`

These files are intended to be opened in ParaView.

## Notes on solver choice

- `red_black_gauss_seidel` is the simplest default and works well for many
  cases.
- `cg` is available on CPU and GPU.
- `miccg0` uses an incomplete Cholesky preconditioner on the CPU when Eigen3 is
  available.
- Particle methods use CFL-based substepping through `max_cfl`, which is useful
  for fast free-fall and impact cases.

## Repository layout

```text
src/
  core/                 core data structures, IO, parameters
  solvers/              SL, PIC, FLIP, APIC solvers
  main.cpp              executable entry point
test/
  SL/                   Semi-Lagrangian sample scenes
  PIC/                  PIC sample scenes
  FLIP/                 FLIP sample scenes
  APIC/                 APIC sample scenes
results/                default output location used by sample configs
PostPro/                optional PICM-PostPro submodule with Slurm/report tools
```

## Post-Processing / Slurm Tools

The report and server scripts live in the separate `PICM-PostPro` repository.
When it is checked out as `PostPro`, use:

```bash
make -C PostPro build
make -C PostPro sbatch
make -C PostPro postpro
make -C PostPro plot
make -C PostPro clean
```

Post-processing outputs stay inside `PostPro/`:

- `PostPro/data/`: CSV files used by plotting.
- `PostPro/data/misc/`: generated configs/logs and temporary post-processing
  inputs.
- `PostPro/img/`: generated `png`, `svg`, `pdf`, and `jpg` figures.

`clean` removes build folders, raw simulation files, and regenerated images
while keeping CSV data.
