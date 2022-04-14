cmake_minimum_required(VERSION 3.3)
project(reformulate)

### c++11 standards
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11 -I/usr/local/include")

### NTL ###
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GMP REQUIRED)
find_package(NTL REQUIRED)


### SCIP ###
if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()


set(CMAKE_LIBRARY_PATH ${SCIP_LIBRARIES})
include_directories(${SCIP_INCLUDE_DIRS})


### EXECUTABLE ###
add_executable(reformulate
   src/main.cpp
   src/EventhdlrReformulate.cpp)
target_link_libraries(reformulate ${NTL_LIBRARY} ${GMP_LIBRARIES} ${M_LIB} ${SCIP_LIBRARIES} ${SCIP_LIB})