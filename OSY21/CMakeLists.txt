cmake_minimum_required(VERSION 3.0.0)
project(OSY21 VERSION 0.1.0)

include(CTest)
enable_testing() 

add_executable(OSY21 main.cpp tool.h)
target_link_libraries(OSY21 ntl m gmp pthread)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
