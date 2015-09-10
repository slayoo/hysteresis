set(f "invfft")
add_executable(${f} ${f}.cpp)

find_package(Boost QUIET COMPONENTS system iostreams)
target_link_libraries(${f} ${Boost_LIBRARIES})
target_include_directories(${f} PUBLIC ${BOOST_INCLUDE_DIRS})

find_package(GSL QUIET)
target_link_libraries(${f} ${GSL_LIBRARIES})
target_include_directories(${f} PUBLIC ${GSL_INCLUDE_DIRS})

add_test(${f} ${f})
