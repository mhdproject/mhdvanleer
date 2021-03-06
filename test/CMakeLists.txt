cmake_minimum_required(VERSION 2.8.7)
set(PROJECT_NAME_STR mhdvanleer)
project(${PROJECT_NAME_STR} C CXX)
find_package(HDF5) # Find non-cmake built HDF5
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
include_directories(/opt/local/include/)

SET(H5FLAGS "-DH5_USE_18_API")
add_definitions(${H5FLAGS})

include_directories(${HDF5_INCLUDE_PATH})


find_package(Threads REQUIRED)

if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-Wall -ansi -Wno-deprecated -pthread)
endif()

if(MSVC)
    #vc 2012 fix for vararg templates
    set(MSVC_COMPILER_DEFS "-D_VARIADIC_MAX=10")
endif() 
#-------------------
# set common include folder for module
#-------------------
set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)

#-------------------
# Module source
#-------------------
set(PROJECT_LIB_NAME ${PROJECT_NAME_STR}_impl)

include_directories(${COMMON_INCLUDES})
set(SRC_FILES
        ${PROJECT_SOURCE_DIR}/src/boundary.cpp ${PROJECT_SOURCE_DIR}/src/cooling.cpp ${PROJECT_SOURCE_DIR}/src/eigenvectors.cpp ${PROJECT_SOURCE_DIR}/src/flux.cpp ${PROJECT_SOURCE_DIR}/src/hlld.cpp ${PROJECT_SOURCE_DIR}/src/initialise.cpp ${PROJECT_SOURCE_DIR}/src/initialise_jet.cpp ${PROJECT_SOURCE_DIR}/src/locate.cpp ${PROJECT_SOURCE_DIR}/src/log.cpp ${PROJECT_SOURCE_DIR}/src/maxspeed.cpp ${PROJECT_SOURCE_DIR}/src/molcool.cpp ${PROJECT_SOURCE_DIR}/src/output.cpp ${PROJECT_SOURCE_DIR}/src/riemann.cpp ${PROJECT_SOURCE_DIR}/src/roe.cpp ${PROJECT_SOURCE_DIR}/src/sgn.cpp ${PROJECT_SOURCE_DIR}/src/tabfind.cpp ${PROJECT_SOURCE_DIR}/src/update.cpp ${PROJECT_SOURCE_DIR}/src/vanleer.cpp)
add_library(${PROJECT_LIB_NAME} ${SRC_FILES})

add_subdirectory(${EXT_PROJECTS_DIR}/gtest)

#-------------------
# Test
#-------------------
add_definitions(${MSVC_COMPILER_DEFS})
enable_testing()
set(PROJECT_TEST_NAME ${PROJECT_NAME_STR}_test)
include_directories(${GTEST_INCLUDE_DIRS} ${COMMON_INCLUDES})

file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)
add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES})
add_dependencies(${PROJECT_TEST_NAME} googletest)
add_executable(main ${PROJECT_SOURCE_DIR}/src/main.cpp ${SRC_FILES})
if (APPLE)
    target_link_libraries(${PROJECT_TEST_NAME}
            ${PROJECT_LIB_NAME}
            ${MPI_CXX_LIBRARIES}
            /opt/local/lib/libhdf5.dylib
            /opt/local/lib/libz.dylib
            ${GTEST_LIBS_DIR}/libgtest.a
            ${GTEST_LIBS_DIR}/libgtest_main.a
            )
    target_link_libraries(main
            ${PROJECT_LIB_NAME}
            ${MPI_CXX_LIBRARIES}
            /opt/local/lib/libhdf5.dylib
            /opt/local/lib/libz.dylib
            )
endif ()

if (NOT WIN32 OR MINGW OR APPLE)
    target_link_libraries(${PROJECT_TEST_NAME}
		${PROJECT_LIB_NAME}
            ${MPI_CXX_LIBRARIES}
            ${HDF5_CXX_LIBRARIES}
            ${HDF5_C_LIBRARIES}
            ${HDF5_LIBRARIES}
            #  /opt/local/lib/libhdf5.dylib
            # /opt/local/lib/libz.dylib
        ${GTEST_LIBS_DIR}/libgtest.a
        ${GTEST_LIBS_DIR}/libgtest_main.a
        )
    target_link_libraries(main
            ${PROJECT_LIB_NAME}
            ${MPI_CXX_LIBRARIES}
            ${HDF5_CXX_LIBRARIES}
            ${HDF5_C_LIBRARIES}
            ${HDF5_LIBRARIES}
            #  /opt/local/lib/libhdf5.dylib
            # /opt/local/lib/libz.dylib
            )
else()
    target_link_libraries(${PROJECT_TEST_NAME}
        debug ${GTEST_LIBS_DIR}/DebugLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest${CMAKE_FIND_LIBRARY_SUFFIXES}
        optimized ${GTEST_LIBS_DIR}/ReleaseLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest${CMAKE_FIND_LIBRARY_SUFFIXES}
        )
    target_link_libraries(${PROJECT_TEST_NAME}
        debug ${GTEST_LIBS_DIR}/DebugLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main${CMAKE_FIND_LIBRARY_SUFFIXES}
        optimized ${GTEST_LIBS_DIR}/ReleaseLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main${CMAKE_FIND_LIBRARY_SUFFIXES}
        )
endif()

target_link_libraries(${PROJECT_TEST_NAME} ${CMAKE_THREAD_LIBS_INIT})

add_test(test1 ${PROJECT_TEST_NAME})
