set(f "rhloop")

add_test(NAME ${f}_init COMMAND "ncgen" "-k" "hdf5" "${CMAKE_CURRENT_SOURCE_DIR}/${f}.cdl")
add_test(NAME ${f}_calc COMMAND "${CMAKE_BINARY_DIR}/main" "${CMAKE_CURRENT_BINARY_DIR}/${f}.nc")

add_executable(${f}_plot ${f}_plot.cpp)
target_link_libraries(${f}_plot ${Boost_LIBRARIES};${HDF5_LIBRARIES})
target_include_directories(${f}_plot PUBLIC ${BOOST_INCLUDE_DIRS};${HDF5_INCLUDE_DIRS})
add_test(${f}_plot ${f}_plot)
