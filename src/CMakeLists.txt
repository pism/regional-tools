cmake_minimum_required (VERSION 3.1)
project(BasinCalculator CXX)

find_package(GSL REQUIRED)

add_library(basincalculator
  DEM.cc
  accumulated_flow.cc
  initialize_mask.cc
  upslope_area.cc)

target_include_directories(basincalculator PUBLIC ${GSL_INCLUDE_DIRS})
target_link_libraries(basincalculator ${GSL_LIBRARY})
