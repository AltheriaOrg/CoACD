# Logia plugin build

This fork tracks upstream [SarahWeiii/CoACD](https://github.com/SarahWeiii/CoACD)
unchanged. The only addition is this folder, which builds `coacd_lib.dll`: the
Unity plugin Logia loads from `Assets/Various/Plugins/CoACD/Release/`.

## Why a shim instead of upstream's `lib_coacd.dll`

Logia's P/Invoke bindings (`MeshProcess.CoACD`, `DebugCPP`) were written
against the 2022 fork, which exported its own seven functions. `coacd_lib.cpp`
re-exports exactly those names and signatures on top of upstream's
`coacd::CoACD()`, so a new DLL drops into Logia with no C# change:

| Export | Does |
|---|---|
| `LoadMesh(float* xyz, uint nPts, uint* tris, uint nTris)` | copies a Unity mesh into a `coacd::Mesh` |
| `DestroyMesh(mesh)` | frees it |
| `GenerateConvexDecomposition(mesh, threshold)` | runs CoACD; never returns null, never throws |
| `GetNConvexHulls(hulls)` / `GetConvexHull(hulls, i, out)` | read the result |
| `DestroyConvexHulls(hulls)` | frees the result and every buffer its hulls point into |
| `RegisterDebugCallback(cb)` | routes CoACD's log (warnings and errors) to Unity's console |

The `ConvexHull` struct matches Logia's C# declaration field for field.

## Parameters

Logia passes only `threshold` (the collider quality: 0.5, 0.2, 0.1, 0.05 or
0.02). Everything else is fixed in `coacd_lib.cpp`:

- `threshold` clamped to 0.01–1, as the old fork did.
- `preprocess = "auto"`: remesh with OpenVDB only when the input is not
  manifold. The old fork always ran its manifold step.
- `pca = true`, as the old fork did (upstream's default is off).
- `seed = 0`: the same model gives the same colliders on every import. The old
  fork seeded from the clock.
- OpenMP uses all cores but one, leaving one for the Editor's main thread.
- Everything else is upstream's default.

Changing any of these changes the colliders new imports get. Existing
projects are unaffected: their hulls are stored in the `.lmd`.

## Build (Windows)

Visual Studio 2022 or later with the C++ workload. The first configure fetches
Boost, OpenVDB, zlib and spdlog.

```
cmake -S . -B build -G "Visual Studio 18 2026" -A x64 -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded -DOPENVDB_CORE_SHARED=OFF -DTBB_TEST=OFF -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build --config Release --target coacd_lib logia_smoke_test
```

The flags are upstream's own Windows CI flags. `cmake` must be on `PATH`
(OpenVDB's patch step calls it by name); from a plain shell, add Visual
Studio's `Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin`.

Output: `build/logia/Release/coacd_lib.dll`. The C++ runtime is linked
statically; the only runtime DLL it needs is `VCOMP140.DLL` (OpenMP), which
the old build needed too.

## Check before copying into Logia

```
cd build/logia/Release
logia_smoke_test.exe 8 0.05
```

It loads the DLL by name, decomposes a 14,400-triangle torus eight times and
prints hulls and private memory per run. Hull counts must be identical and
memory flat from run 2 (it settles around 20–22 MB, with early runs spiking to ~45 MB as the
allocator keeps and returns pages, with no upward trend).

Reference numbers on this test (threshold 0.05, same machine):

| Build | Per run | Hulls | Memory |
|---|---|---|---|
| 2022 fork (Logia before this) | ~21 s | varies per run (clock seed) | grows, 82.6 MB after two runs |
| This build | ~12.5 s | 10 every run | flat |

A build that takes ~85 s per run was compiled without optimisation: check
the MSVC branch of `CMAKE_CXX_FLAGS_RELEASE` in the root `CMakeLists.txt`.

## Syncing with upstream again

```
git fetch upstream
git merge upstream/main
```

Upstream never touches `logia/`. Outside it, the fork changes two things in
the root `CMakeLists.txt`: the `add_subdirectory(logia)` block at the end,
and MSVC Release flags (`/O2 /Ob2 /DNDEBUG` instead of `-O3`, which MSVC
ignores). The second is an upstream bug worth sending upstream; once it is
fixed there, drop ours in the merge.
