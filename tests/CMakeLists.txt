set(CMAKE_CXX_FLAGS "-I${CMAKE_SOURCE_DIR}/src ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

add_subdirectory(invfft)
