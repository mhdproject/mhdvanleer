cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(PROJECT_NAME_STR mhdvanleer)
project(${PROJECT_NAME_STR} C CXX)
find_package(HDF5) # Find non-cmake built HDF5
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
include_directories(/opt/local/include/)
include_directories(${PROJECT_SOURCE_DIR}/include/)
include_directories(.)


SET(H5FLAGS "-DH5_USE_18_API")
add_definitions(${H5FLAGS})

include_directories(${HDF5_INCLUDE_PATH})


find_package(Threads REQUIRED)

if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-Wall -Wextra -ansi -Wno-deprecated -pthread -pedantic)
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
        ${PROJECT_SOURCE_DIR}/include/TwoDim.h
        ${PROJECT_SOURCE_DIR}/include/boundary.h
        ${PROJECT_SOURCE_DIR}/include/cooling.h
        ${PROJECT_SOURCE_DIR}/include/eigenvectors.h
        ${PROJECT_SOURCE_DIR}/include/flux.h
        ${PROJECT_SOURCE_DIR}/include/global.h
        ${PROJECT_SOURCE_DIR}/include/hlld.h
        ${PROJECT_SOURCE_DIR}/include/initialise.h
        ${PROJECT_SOURCE_DIR}/include/locate.h
        ${PROJECT_SOURCE_DIR}/include/log.h
        ${PROJECT_SOURCE_DIR}/include/main.h
        ${PROJECT_SOURCE_DIR}/include/maxspeed.h
        ${PROJECT_SOURCE_DIR}/include/molcool.h
        ${PROJECT_SOURCE_DIR}/include/output.h
        ${PROJECT_SOURCE_DIR}/include/riemann.h
        ${PROJECT_SOURCE_DIR}/include/roe.h
        ${PROJECT_SOURCE_DIR}/include/sgn.h
        ${PROJECT_SOURCE_DIR}/include/tabfind.h
        ${PROJECT_SOURCE_DIR}/include/tnt.h
        ${PROJECT_SOURCE_DIR}/include/tnt_array1d.h
        ${PROJECT_SOURCE_DIR}/include/tnt_array1d_utils.h
        ${PROJECT_SOURCE_DIR}/include/tnt_array2d.h
        ${PROJECT_SOURCE_DIR}/include/tnt_array2d_utils.h
        ${PROJECT_SOURCE_DIR}/include/tnt_array3d.h
        ${PROJECT_SOURCE_DIR}/include/tnt_array3d_utils.h
        ${PROJECT_SOURCE_DIR}/include/tnt_cmat.h
        ${PROJECT_SOURCE_DIR}/include/tnt_fortran_array2d.h
        ${PROJECT_SOURCE_DIR}/include/tnt_fortran_array2d_utils.h
        ${PROJECT_SOURCE_DIR}/include/tnt_fortran_array3d.h
        ${PROJECT_SOURCE_DIR}/include/tnt_fortran_array3d_utils.h
        ${PROJECT_SOURCE_DIR}/include/tnt_i_refvec.h
        ${PROJECT_SOURCE_DIR}/include/tnt_math_utils.h
        ${PROJECT_SOURCE_DIR}/include/tnt_stopwatch.h
        ${PROJECT_SOURCE_DIR}/include/tnt_subscript.h
        ${PROJECT_SOURCE_DIR}/include/tnt_vec.h
        ${PROJECT_SOURCE_DIR}/include/tnt_version.h
        ${PROJECT_SOURCE_DIR}/include/update.h
        ${PROJECT_SOURCE_DIR}/include/vanleer.h
        ${PROJECT_SOURCE_DIR}/include/vanleer_fvsplit.h
        ${PROJECT_SOURCE_DIR}/src/boundary.cpp
        ${PROJECT_SOURCE_DIR}/src/cooling.cpp
        ${PROJECT_SOURCE_DIR}/src/eigenvectors.cpp
        ${PROJECT_SOURCE_DIR}/src/flux.cpp
        ${PROJECT_SOURCE_DIR}/src/hlld.cpp
        ${PROJECT_SOURCE_DIR}/src/initialise.cpp
        ${PROJECT_SOURCE_DIR}/src/initialise_jet.cpp
        ${PROJECT_SOURCE_DIR}/src/locate.cpp
        ${PROJECT_SOURCE_DIR}/src/log.cpp
        ${PROJECT_SOURCE_DIR}/src/main.cpp
        ${PROJECT_SOURCE_DIR}/src/maxspeed.cpp
        ${PROJECT_SOURCE_DIR}/src/molcool.cpp
        ${PROJECT_SOURCE_DIR}/src/output.cpp
        ${PROJECT_SOURCE_DIR}/src/riemann.cpp
        ${PROJECT_SOURCE_DIR}/src/roe.cpp
        ${PROJECT_SOURCE_DIR}/src/sgn.cpp
        ${PROJECT_SOURCE_DIR}/src/tabfind.cpp
        ${PROJECT_SOURCE_DIR}/src/update.cpp
        ${PROJECT_SOURCE_DIR}/src/vanleer.cpp
        )
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
target_include_directories(main PUBLIC ${PROJECT_SOURCE_DIR}/include/)
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
