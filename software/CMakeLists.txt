PROJECT(MH)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE (Armadillo REQUIRED ARMA_DONT_USE_WRAPPER)

file (GLOB_RECURSE MH_SOURCES "src/*.c*")
file (GLOB_RECURSE MH_HEADERS "include/*.h*")

set (MH_INCLUDE_DIRS "")

foreach(_headersFile ${MH_HEADERS})
    get_filename_component(_dir ${_headersFile} PATH)
    list (APPEND MH_INCLUDE_DIRS ${_dir})
endforeach()

list(REMOVE_DUPLICATES MH_INCLUDE_DIRS)

include_directories(${MH_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS})

add_executable(MH ${MH_SOURCES}) 
target_link_libraries(MH ${ARMADILLO_LIBRARIES})