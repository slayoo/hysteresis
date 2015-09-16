cmake_minimum_required(VERSION 3.0)
project(hysteresis CXX)
enable_testing() 

find_package(Boost QUIET COMPONENTS system iostreams) 
# TODO: what if not found
find_package(GSL QUIET)
# TODO: what if not found
find_package(HDF5 QUIET COMPONENTS CXX HL)
# TODO: what if not found

add_subdirectory(tests)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wextra")

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "")
else()
  set(CMAKE_CXX_FLAGS_DEBUG "-DBZ_DEBUG -lblitz -Og")
endif()

# TODO: find sundials
# TODO: set sundials_LIBRARIES properly
set(sundials_LIBRARIES "sundials_cvode;sundials_nvecserial")

#find blitz
#find boost: units

set(f "main")
add_executable(${f} ${f}.cpp)
target_link_libraries(${f} ${sundials_LIBRARIES};${GSL_LIBRARIES};${HDF5_LIBRARIES})
unset(f)
