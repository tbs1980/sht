# Similar to the CMake Files in https://github.com/Crascit/DownloadProject
project(sht-test)

include(${CMAKE_SOURCE_DIR}/cmake/DownloadProject.cmake)

download_project(
  PROJ googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG master
  UPDATE_DISCONNECTED 1
)

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})

add_executable(test_sht test_sht.cpp
  test_legendre.cpp
)
target_link_libraries(test_sht gtest gmock_main)

add_custom_target(shttest COMMAND test_sht)