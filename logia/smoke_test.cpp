// Loads coacd_lib.dll the way Unity does (by name, exports only), decomposes
// a torus N times and prints hull counts and process memory after each run.
// Memory must plateau after the first run: that is the leak check.
//
//   logia_smoke_test.exe [runs=8] [threshold=0.05]

#include <windows.h>
#include <psapi.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

struct ConvexHull {
  double *points;
  uint32_t *triangles;
  uint32_t nPoints;
  uint32_t nTriangles;
  double volume;
  double center[3];
  uint32_t meshId;
  double bmin[3];
  double bmax[3];
};

typedef void *(*LoadMeshFn)(const float *, uint32_t, const uint32_t *, uint32_t);
typedef void (*DestroyFn)(void *);
typedef void *(*DecomposeFn)(void *, double);
typedef uint32_t (*CountFn)(void *);
typedef void (*GetHullFn)(void *, uint32_t, ConvexHull *);
typedef void (*Callback)(const char *, int, int);
typedef void (*RegisterFn)(Callback);

static void OnLog(const char *msg, int, int size) {
  std::printf("  log: %.*s\n", size, msg);
}

static double PrivateMB() {
  PROCESS_MEMORY_COUNTERS_EX pmc{};
  GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS *)&pmc,
                       sizeof(pmc));
  return pmc.PrivateUsage / (1024.0 * 1024.0);
}

int main(int argc, char **argv) {
  int runs = argc > 1 ? std::atoi(argv[1]) : 8;
  double threshold = argc > 2 ? std::atof(argv[2]) : 0.05;

  HMODULE dll = LoadLibraryA("coacd_lib.dll");
  if (!dll) {
    std::printf("cannot load coacd_lib.dll (%lu)\n", GetLastError());
    return 1;
  }
  auto loadMesh = (LoadMeshFn)GetProcAddress(dll, "LoadMesh");
  auto destroyMesh = (DestroyFn)GetProcAddress(dll, "DestroyMesh");
  auto decompose = (DecomposeFn)GetProcAddress(dll, "GenerateConvexDecomposition");
  auto count = (CountFn)GetProcAddress(dll, "GetNConvexHulls");
  auto getHull = (GetHullFn)GetProcAddress(dll, "GetConvexHull");
  auto destroyHulls = (DestroyFn)GetProcAddress(dll, "DestroyConvexHulls");
  auto registerCb = (RegisterFn)GetProcAddress(dll, "RegisterDebugCallback");
  if (!loadMesh || !destroyMesh || !decompose || !count || !getHull ||
      !destroyHulls || !registerCb) {
    std::printf("missing export\n");
    return 1;
  }
  registerCb(OnLog);

  // Torus, 120 x 60 segments = 14,400 triangles. Not convex, so CoACD has
  // real work to do.
  const int U = 120, V = 60;
  const float R = 1.0f, r = 0.35f;
  std::vector<float> pts;
  std::vector<uint32_t> tris;
  for (int i = 0; i < U; ++i)
    for (int j = 0; j < V; ++j) {
      float u = 2 * 3.14159265f * i / U, v = 2 * 3.14159265f * j / V;
      pts.push_back((R + r * std::cos(v)) * std::cos(u));
      pts.push_back(r * std::sin(v));
      pts.push_back((R + r * std::cos(v)) * std::sin(u));
    }
  for (int i = 0; i < U; ++i)
    for (int j = 0; j < V; ++j) {
      uint32_t a = i * V + j, b = ((i + 1) % U) * V + j;
      uint32_t c = ((i + 1) % U) * V + (j + 1) % V, d = i * V + (j + 1) % V;
      tris.insert(tris.end(), {a, b, c, a, c, d});
    }

  std::printf("start: %.1f MB private\n", PrivateMB());
  uint32_t first = 0;
  for (int run = 1; run <= runs; ++run) {
    DWORD t0 = GetTickCount();
    void *mesh = loadMesh(pts.data(), (uint32_t)(pts.size() / 3), tris.data(),
                          (uint32_t)(tris.size() / 3));
    void *hulls = decompose(mesh, threshold);
    destroyMesh(mesh);

    uint32_t n = count(hulls), verts = 0;
    for (uint32_t k = 0; k < n; ++k) {
      ConvexHull ch;
      getHull(hulls, k, &ch);
      verts += ch.nPoints;
      for (uint32_t t = 0; t < ch.nTriangles * 3; ++t)
        if (ch.triangles[t] >= ch.nPoints) {
          std::printf("hull %u: index out of range\n", k);
          return 1;
        }
    }
    destroyHulls(hulls);
    if (run == 1)
      first = n;
    std::printf("run %d: %u hulls, %u verts, %lu ms, %.1f MB private\n", run, n,
                verts, GetTickCount() - t0, PrivateMB());
    if (n == 0 || n != first) {
      std::printf("unexpected hull count\n");
      return 1;
    }
  }
  return 0;
}
