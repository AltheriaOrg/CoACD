// Logia's Unity-facing exports, kept ABI-compatible with the 2022 fork's
// coacd_lib.dll so Logia's P/Invoke bindings (MeshProcess.CoACD, DebugCPP)
// work unchanged. Everything here is a thin shim over coacd::CoACD().
//
// Ownership: LoadMesh / DestroyMesh pair up, GenerateConvexDecomposition /
// DestroyConvexHulls pair up. DestroyConvexHulls frees everything the hull
// pointers point into. Nothing else is allocated across the boundary.

#include "coacd.h"
#include "../src/logger.h"

#include <spdlog/sinks/base_sink.h>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <exception>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(_WIN32)
#define LOGIA_API extern "C" __declspec(dllexport)
#else
#define LOGIA_API extern "C" __attribute__((visibility("default")))
#endif

namespace {

// Must match MeshProcess.CoACD.ConvexHull in Logia field for field
// (sequential layout, default packing).
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

struct HullSet {
  std::vector<std::vector<double>> points;
  std::vector<std::vector<uint32_t>> triangles;
  std::vector<ConvexHull> hulls;
};

// ---- logging -------------------------------------------------------------

// Colour indices understood by Logia's DebugCPP.OnDebugCallback.
enum class Color : int { Red, Green, Blue, Black, White, Yellow, Orange };

typedef void (*FuncCallBack)(const char *message, int color, int size);
std::atomic<FuncCallBack> g_callback{nullptr};

void Send(const std::string &msg, Color color) {
  FuncCallBack cb = g_callback.load();
  if (cb)
    cb(msg.c_str(), static_cast<int>(color), static_cast<int>(msg.size()));
}

class UnitySink : public spdlog::sinks::base_sink<std::mutex> {
protected:
  void sink_it_(const spdlog::details::log_msg &msg) override {
    std::string text = "[CoACD] ";
    text.append(msg.payload.data(), msg.payload.size());
    Color color = msg.level >= spdlog::level::err    ? Color::Red
                  : msg.level == spdlog::level::warn ? Color::Yellow
                                                     : Color::White;
    Send(text, color);
  }
  void flush_() override {}
};

void AttachUnitySink() {
  static std::once_flag once;
  std::call_once(once, [] {
    auto logger = coacd::logger::get();
    // No console in the Unity Editor: route everything to the callback.
    logger->sinks().clear();
    logger->sinks().push_back(std::make_shared<UnitySink>());
    // CoACD logs every parameter and progress step at info; warnings are
    // what Logia wants in its console.
    logger->set_level(spdlog::level::warn);
  });
}

} // namespace

// ---- exports ---------------------------------------------------------------

LOGIA_API void RegisterDebugCallback(FuncCallBack cb) {
  g_callback.store(cb);
  AttachUnitySink();
}

LOGIA_API void *LoadMesh(const float *points, uint32_t countPoints,
                         const uint32_t *triangles, uint32_t countTriangles) {
  auto *mesh = new coacd::Mesh();
  mesh->vertices.resize(countPoints);
  for (uint32_t i = 0; i < countPoints; ++i)
    mesh->vertices[i] = {points[i * 3], points[i * 3 + 1], points[i * 3 + 2]};
  mesh->indices.resize(countTriangles);
  for (uint32_t i = 0; i < countTriangles; ++i)
    mesh->indices[i] = {(int)triangles[i * 3], (int)triangles[i * 3 + 1],
                        (int)triangles[i * 3 + 2]};
  return mesh;
}

LOGIA_API void DestroyMesh(void *mesh) {
  delete static_cast<coacd::Mesh *>(mesh);
}

// Never returns null: on failure the error goes to the debug callback and
// the result is an empty set, so Logia's GetNConvexHulls/DestroyConvexHulls
// stay valid. No C++ exception may cross into Mono.
LOGIA_API void *GenerateConvexDecomposition(void *mesh, double threshold) {
  auto *result = new HullSet();
  try {
    AttachUnitySink();

#ifdef _OPENMP
    // Leave one core for the Editor's main thread, as the old fork did.
    int cores = (int)std::thread::hardware_concurrency();
    omp_set_num_threads(std::max(1, cores - 1));
#endif

    const auto &input = *static_cast<coacd::Mesh *>(mesh);
    std::vector<coacd::Mesh> parts = coacd::CoACD(
        input,
        std::clamp(threshold, 0.01, 1.0), // same clamp as the old fork
        /*max_convex_hull*/ -1,
        /*preprocess*/ "auto", // remesh only when the input isn't manifold
        /*prep_resolution*/ 50,
        /*sample_resolution*/ 2000,
        /*mcts_nodes*/ 20,
        /*mcts_iteration*/ 150,
        /*mcts_max_depth*/ 3,
        /*pca*/ true, // the old fork always ran PCA
        /*merge*/ true,
        /*decimate*/ false,
        /*max_ch_vertex*/ 256,
        /*extrude*/ false,
        /*extrude_margin*/ 0.01,
        /*apx_mode*/ "ch",
        /*seed*/ 0, // fixed: the same model gives the same colliders
        /*real_metric*/ false);

    const size_t n = parts.size();
    result->points.resize(n);
    result->triangles.resize(n);
    result->hulls.resize(n);
    for (size_t i = 0; i < n; ++i) {
      const coacd::Mesh &part = parts[i];
      auto &pts = result->points[i];
      auto &tris = result->triangles[i];
      pts.reserve(part.vertices.size() * 3);
      for (const auto &v : part.vertices)
        pts.insert(pts.end(), {v[0], v[1], v[2]});
      tris.reserve(part.indices.size() * 3);
      for (const auto &t : part.indices)
        tris.insert(tris.end(),
                    {(uint32_t)t[0], (uint32_t)t[1], (uint32_t)t[2]});

      ConvexHull ch{};
      ch.points = pts.data();
      ch.triangles = tris.data();
      ch.nPoints = (uint32_t)part.vertices.size();
      ch.nTriangles = (uint32_t)part.indices.size();
      ch.meshId = (uint32_t)i;
      result->hulls[i] = ch;
    }
  } catch (const std::exception &e) {
    Send(std::string("[CoACD] decomposition failed: ") + e.what(), Color::Red);
    result->hulls.clear();
  } catch (...) {
    Send("[CoACD] decomposition failed: unknown error", Color::Red);
    result->hulls.clear();
  }
  return result;
}

LOGIA_API uint32_t GetNConvexHulls(void *hulls) {
  return hulls ? (uint32_t)static_cast<HullSet *>(hulls)->hulls.size() : 0;
}

LOGIA_API void GetConvexHull(void *hulls, uint32_t index, ConvexHull *out) {
  *out = static_cast<HullSet *>(hulls)->hulls[index];
}

LOGIA_API void DestroyConvexHulls(void *hulls) {
  delete static_cast<HullSet *>(hulls);
}
