# CMakeList.txt : CMake project for NbodyTreeCode, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)

# Add source to this project's executable.
add_executable (NbodyTreeCode "NbodyTreeCode.cpp" "NbodyTreeCode.h" "Body.h" "Rect.h" "AdaptiveQuadtree.h" "AdaptiveQuadtree.cpp")

# TODO: Add tests and install targets if needed.
