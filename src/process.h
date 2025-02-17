#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include <math.h>
#include <limits>
#include <typeinfo>
#include <algorithm>
#include <assert.h>
#include <regex>

#include "io.h"
#include "clip.h"
#include "config.h"
#include "model_obj.h"
#include "cost.h"

#pragma once
void ManifoldPreprocess(Params& params, Model& input, Model& output);
void ManifoldPreprocess(Params &params, Model &m, ofstream &of);
void MergeCH(Model &ch1, Model &ch2, Model &ch);
double MergeConvexHulls(Model &m, vector<Model> &meshs, vector<Model> &cvxs, Params &params, ofstream &of, double epsilon = 0.02, double threshold = 0.01);
void Compute(ofstream &of, Model &mesh, Params &params);
vector<Model>* ComputeLib(ofstream& of, Model& mesh, Params& params);

inline void addNeighbor(map<pair<int, int>, pair<int, int>> &edge_map, pair<int, int> &edge, vector<int> &neighbors, int idx)
{
  int first = edge_map[edge].first;
  int second = edge_map[edge].second;
  if (first != idx && first != -1)
    neighbors.push_back(first);
  if (second != idx && second != -1)
    neighbors.push_back(second);
}

inline int32_t FindMinimumElement(const vector<double> d, double *const m, const int32_t begin, const int32_t end)
{
  int32_t idx = -1;
  double min = (std::numeric_limits<double>::max)();
  for (size_t i = begin; i < size_t(end); ++i)
  {
    if (d[i] < min)
    {
      idx = i;
      min = d[i];
    }
  }

  *m = min;
  return idx;
}
