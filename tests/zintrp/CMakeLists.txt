set(f "zintrp")
add_executable(${f} ${f}.cpp)
target_link_libraries(${f} ${Boost_LIBRARIES};${GSL_LIBRARIES})
target_include_directories(${f} PUBLIC ${BOOST_INCLUDE_DIRS};${GSL_INCLUDE_DIRS})
add_test(${f} ${f})
