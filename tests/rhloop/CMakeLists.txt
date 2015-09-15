set(f "rhloop")

add_executable(${f}_plot ${f}_plot.cpp)
target_link_libraries(${f}_plot ${Boost_LIBRARIES};${HDF5_LIBRARIES})
target_include_directories(${f}_plot PUBLIC ${BOOST_INCLUDE_DIRS};${HDF5_INCLUDE_DIRS})
add_test(${f}_plot ${f}_plot)
