project(bench-sht CXX)
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
add_executable(sht_gsl sht_gsl.cpp)
target_link_libraries(sht_gsl ${GSL_LIBRARIES})

add_executable(sht_multiprecision sht_multiprecision.cpp)