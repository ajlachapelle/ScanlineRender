cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(ScanlineRender)
SET(VTK_DIR ./VTK-8.1.2)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(rasterizer rasterizer)

if(VTK_LIBRARIES)
  target_link_libraries(rasterizer ${VTK_LIBRARIES})
else()
  target_link_libraries(rasterizer vtkHybrid)
endif()

