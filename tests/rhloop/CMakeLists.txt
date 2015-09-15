set(f "rhloop")

add_executable(${f} ${f}.cpp)
target_link_libraries(${f} ${Boost_LIBRARIES};${HDF5_LIBRARIES})
target_include_directories(${f} PUBLIC ${BOOST_INCLUDE_DIRS};${HDF5_INCLUDE_DIRS})
add_test(${f} ${f})
