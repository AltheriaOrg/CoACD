/*
Copyright (c) 2018 Jingwei Huang. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to the authors of this software, without
imposing a separate written license agreement for such Enhancements, then you
hereby grant the following license: a non-exclusive, royalty-free perpetual
license to install, use, modify, prepare derivative works, incorporate into
other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.
*/
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <vector>

#include "model.h"
#include "../model_obj.h"

extern int g_sharp;
int Manifold(ofstream& of, string input_model, Model& output, int resolution)
{
  clock_t start, end;

  Model_OBJ obj;
  char *cstr = new char[input_model.length() + 1];
  strcpy(cstr, input_model.c_str());
  obj.Load(cstr);
  delete [] cstr;
  
  of << " - Pre-processing (manifold)" << endl;
  of << "\tresolution: " << resolution << endl;
  cout << " - Pre-processing (manifold)" << endl;
  cout << "\tresolution: " << resolution << endl;
  
  start = clock();
  obj.Process_Manifold(resolution);
  end = clock();
  of << "Manifold time: " << double(end-start)/CLOCKS_PER_SEC << "s" << endl;
  cout << "Manifold time: " << double(end-start)/CLOCKS_PER_SEC << "s" << endl;

  output.Load(obj.vertices, obj.face_indices);
  
  return 0; 
}

int Manifold(Model& input, Model& output, int resolution)
{
	clock_t start, end;
	Model_OBJ obj;
	for (int i = 0; i < input.points.size(); ++i)
		obj.vertices.push_back(glm::dvec3(input.points[i][0],input.points[i][1] ,input.points[i][2]));
	for (int i = 0; i < input.triangles.size(); ++i)
		obj.face_indices.push_back(glm::ivec3(input.triangles[i][0], input.triangles[i][1] , input.triangles[i][2]));
	obj.Process_Manifold(resolution);

	output.Load(obj.vertices, obj.face_indices);

	return 0;
}
