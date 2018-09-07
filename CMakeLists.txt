cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(polynomial_curve_fit)
# find OpenCV
FIND_PACKAGE( OpenCV REQUIRED )
  if(OpenCV_FOUND)
    MESSAGE(STATUS "Found OpenCV libraries at ${OpenCV_LIBRARIES} and include files at ${OpenCV_INCLUDE_DIRS}")
	include_directories(SYSTEM ${OpenCV2_INCLUDE_DIRS})
	link_directories(${OpenCV2_LIBRARY_DIRS})
	link_directories(C:/OpenCV/x64/vc12/lib)
  endif(OpenCV_FOUND)


add_executable (polynomial_curve_fit 
  curve.h
  curve.cpp
  Test_curve.cpp
)
#MRF
target_link_libraries (polynomial_curve_fit ${OpenCV_LIBS})