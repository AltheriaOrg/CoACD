#include <algorithm>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <assert.h>
#include <algorithm>
#include <set>
#include <map>
#include <unordered_map>

#pragma once

#include "shape.h"
#include "model_obj.h"

//////////////// IO ////////////////
void SaveConfig(ofstream &of, Params params);
void SaveOBJ(const string &filename, vector<Model> parts, Params &params);
void SaveOBJwoRecover(const string& filename, vector<Model> parts);
void SaveOBJS(const string foldername, const string &filename, vector<Model> parts, Params &params);
bool WriteVRML(ofstream &fout, Model mesh);
void SaveVRML(const string &fileName, vector<Model> meshs, Params &params);
