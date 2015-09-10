cmake_minimum_required(VERSION 3.0)
project(hysteresis CXX)
enable_testing() 
add_subdirectory(tests)




add_executable(hysteresis test.cpp)

# TODO: Ofast etc
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

# TODO: find sundials
# TODO: set sundials_LIBRARIES properly
set(sundials_LIBRARIES "sundials_cvode;sundials_nvecserial")

#find blitz
#find boost: units

target_link_libraries(hysteresis ${sundials_LIBRARIES})

