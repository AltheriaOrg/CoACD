# Logia's Unity plugin: coacd_lib.dll with the exports Logia's P/Invoke
# bindings expect. See logia/README.md.
add_library(coacd_lib SHARED coacd_lib.cpp ${PROJECT_SOURCE_DIR}/public/coacd.cpp)
target_link_libraries(coacd_lib PRIVATE coacd_objects coacd_interface
                                        openvdb_static spdlog::spdlog)
target_include_directories(coacd_lib PRIVATE ${PROJECT_SOURCE_DIR}/public)
# coacd.h would otherwise mark the C++ API dllimport.
target_compile_definitions(coacd_lib PRIVATE COACD_STATIC)
# Upstream sets a global "lib" prefix on MSVC; Logia loads "coacd_lib".
set_target_properties(coacd_lib PROPERTIES PREFIX "")

# Loads coacd_lib.dll by name like Unity does and checks memory stays flat
# across repeated decompositions.
if(WIN32)
    add_executable(logia_smoke_test smoke_test.cpp)
    target_link_libraries(logia_smoke_test PRIVATE psapi)
    add_dependencies(logia_smoke_test coacd_lib)
endif()
