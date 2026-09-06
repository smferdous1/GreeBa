cmake_minimum_required(VERSION 3.10)
project(GreeBa)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

# Some SCIP packages export libscip without the namespace OR-Tools expects.
# Forward its usage requirements before OR-Tools runs its warning-producing fallback.
find_package(SCIP CONFIG QUIET)
if(SCIP_FOUND AND TARGET libscip AND NOT TARGET SCIP::libscip)
    add_library(SCIP::libscip INTERFACE IMPORTED)
    set_target_properties(SCIP::libscip PROPERTIES INTERFACE_LINK_LIBRARIES libscip)
endif()

# Use the installation supplied by the environment, CMAKE_PREFIX_PATH, or ortools_DIR.
# OR-Tools' package configuration locates its own dependency modules.
find_package(ortools CONFIG REQUIRED)
message(STATUS "OR-Tools CMake package: ${ortools_DIR}")

include(CTest)

add_subdirectory(src)
add_subdirectory(app)
