cmake_minimum_required(VERSION 3.10)
project(GreeBa)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

list(APPEND CMAKE_PREFIX_PATH /opt/homebrew)
list(APPEND CMAKE_MODULE_PATH /opt/homebrew/Cellar/or-tools/9.15_6/lib/cmake/ortools/modules)
find_package(ortools CONFIG REQUIRED)

include(CTest)

add_subdirectory(src)
add_subdirectory(app)
