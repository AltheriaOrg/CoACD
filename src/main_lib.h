#pragma once

#ifdef WIN32
#define EXTERN  extern "C" __declspec(dllexport)
#else
#define EXTERN  extern "C"
#endif

EXTERN void* GenerateConvexDecomposition(Model*, double);