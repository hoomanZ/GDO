# This is the root Finite-Element CMakeLists file.

#----------------------------------------------------------------------------
# Create the project and check cmake version
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(FiniteElement CXX)

set(MAJOR_VERSION "1")
set(MINOR_VERSION "0")
set(PATCH_VERSION "0")

set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_BUILD_TYPE "Release")
set(BUILD_SHARED_LIBS ON)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -Wall -llapack -lblas -larpack -lgfortran")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -Wall -O2 -larmadillo")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -Wall -I usr/lib -O1  -larmadillo -llapack -lblas -lgfortran -lm")
##set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -Wall -I /usr/include/armadilo -O2 -llapack -larpack -lblas")
##set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -Wall -I -O2 -llapack -larpack -lblas")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -Wall -O2 -L/user/include/armadillo-4.450.4/include -llapack -lblas")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -Wall -O2 -L/user/include/armadillo-4.450.4/include -#DARMA_DONT_USE_WRAPPER -llapack -lblas -lgfortran -lm")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -Wall -O2 -L/home/FiniteElement/src/armadillo-4.500.0/include -DARMA_DONT_USE_WRAPPER -lgfortran -lblas -llapack")
message("ENV_FLAGS = $ENV{FFLAGS}")
message("CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

set(BASEPATH ${CMAKE_SOURCE_DIR})
include_directories(${BASEPATH})

include("/usr/share/Armadillo/CMake/ArmadilloLibraryDepends.cmake")

#----------------------------------------------------------------------------
# Set module path to search for local Find<XXX>.cmake files
#----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/Modules/")

#set(ARMADILLO_INCLUDE_DIRS "/usr/include/armadillo-4.450.4/include/")
#set(ARMADILLO_LIBRARY "/usr/share/doc/")

#set(ARMADILLO_INCLUDE_DIRS "/usr/include")
#set(ARMADILLO_LIBRARY "/usr/lib")


#----------------------------------------------------------------------------
# Needs boost libraries: at least 1.49
#----------------------------------------------------------------------------
find_package(Boost 1.49 COMPONENTS REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Needs boost libraries: at least 1.49
#----------------------------------------------------------------------------
find_package(Armadillo  REQUIRED)
include_directories(${Armadillo_INCLUDE_DIRS})


#----------------------------------------------------------------------------
# Test variant
#----------------------------------------------------------------------------
add_executable(finiteElement finiteElement.cc)

target_link_libraries(finiteElement
  FiniteElement_LIB  
  ${Boost_LIBRARIES})
#  ${Armadillo_LIBRARIES})

add_executable(testMPM testMPM.cc)

target_link_libraries(testMPM
  FiniteElement_LIB  
  ${Boost_LIBRARIES})

add_executable(test_peri test_peri.cc)

target_link_libraries(test_peri
  FiniteElement_LIB  
  ${Boost_LIBRARIES})
#  ${Armadillo_LIBRARIES})

add_executable(test_state test_state.cc)

target_link_libraries(test_state
  FiniteElement_LIB  
  ${Boost_LIBRARIES})
#  ${Armadillo_LIBRARIES})

add_executable(test_derivative test_derivative.cc)

target_link_libraries(test_derivative
  FiniteElement_LIB  
  ${Boost_LIBRARIES})
#  ${Armadillo_LIBRARIES})
#----------------------------------------------------------------------------
# Set up library sources
#----------------------------------------------------------------------------
set(FiniteElement_SRC
  Tensor.cc
  Node.cc
  TrilinearVolumeElement.cc
  SymmetricMaterial.cc
  SymmetricMaterialTrilinearElement.cc
  Body.cc
  Box.cc
  ComplicatedGeometry.cc
  Solver.cc
  MaterialPoint.cc
  MPMBox.cc
  MPMBoxPhases.cc
  PeriPMBMaterial.cc
  PeriMaterialPoint.cc
  PeriBond.cc
  PeriBox.cc
  PeriCrack.cc
  Peridynamics.cc
  MPMSolver.cc

  StateMaterialPoint.cc
  StateBond.cc

  derivative_operator.cc
)
set(GEOMETRY_MATH_SRC
 # GeometryMath/Box3D.cc
   GeometryMath/IntVector3D.cc
   GeometryMath/Point3D.cc
   GeometryMath/Polygon3D.cc
   GeometryMath/Vector3D.cc
   GeometryMath/Matrix3D.cc
   GeometryMath/Matrix.cc
 )
#set(GEOMETRY_PIECE_SRC
#  GeometryPiece/BoxGeometryPiece.cc
#  GeometryPiece/GeometryPiece.cc
#  GeometryPiece/GeometryPieceFactory.cc
#  GeometryPiece/GeometryReader.cc
#)
#set(SHAPE_FUNCTION_SRC
#  ShapeFunctions/MPMShapeFunctionFactory.cc
 # ShapeFunctions/LinearShapeFunction.cc
#)
#set(CONTACT_SRC
 # Contact/MPMContactFactory.cc
#  Contact/MPMContact.cc
#  Contact/MPMFreeContact.cc
#  Contact/MPMFrictionContact.cc
#  Contact/MPMFrictionlessContact.cc
#  Contact/MPMVelocityContact.cc
#)

#add_library(BrMPM_LIB 
add_library(FiniteElement_LIB 
   ${FiniteElement_SRC}
   ${GEOMETRY_MATH_SRC})

#----------------------------------------------------------------------------
# Test input file reader
#----------------------------------------------------------------------------
#foreach(dir ${PROBLEMSPEC_INCLUDE_DIR})
# message(STATUS "Dir = ${dir}")
#endforeach()

#add_executable(test_reader test_reader.cc)
#target_link_libraries(test_reader 
#  BrMPM_LIB
#  ${PROBLEMSPEC_LIBRARY}
#  ${Boost_LIBRARIES}
#  ${TRIANGLE_LIBRARY}
#  ${VTK_LIBRARIES}
#  ${XML2_LIBRARY})

#----------------------------------------------------------------------------
# Test mpm
#----------------------------------------------------------------------------
# add_executable(test_two_contact test_two_contact.cc)
# target_link_libraries(test_two_contact 
#   BrMPM_LIB
#   ${PROBLEMSPEC_LIBRARY}
#   ${Boost_LIBRARIES}
#   ${VTK_LIBRARIES}
#   ${XML2_LIBRARY})

