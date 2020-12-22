cmake_minimum_required(VERSION 3.17)
project(TestSofaExecutable)

add_executable(${PROJECT_NAME} main.cpp)

set(CMAKE_CXX_STANDARD ++17)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_PREFIX_PATH "C:/sofa_libs/eigen-3.3.7")

if (WIN32)
    list(APPEND CMAKE_PREFIX_PATH "C:/sofa_libs/eigen-3.3.7")
    list(APPEND CMAKE_PREFIX_PATH "C:/sofa_libs/boost_1_72_0")
    set(BOOST_ROOT "C:/sofa_libs/boost_1_72_0")
    list(APPEND CMAKE_PREFIX_PATH "C:/sofa_libs/sofa")
    list(APPEND CMAKE_PREFIX_PATH "C:/sofa_libs/sofa/lib/win64")

    set(sofa_rooot "C:/Users/Somers/git/sofa-latest/build/install")
    set(sofa_diiir C:/Users/Somers/git/sofa-latest/build/install/include)

    find_package(Boost REQUIRED)
    message("Including Boost header directory:\t ${Boost_INCLUDE_DIRS}/boost_1_72_0" )
    include_directories("${Boost_INCLUDE_DIRS}/boost_1_72_0")
else()
    set(sofa_rooot ~git/sofa/build/install)
    set(sofa_diiir ~/git/sofa/build/install/include)

    find_package(Boost REQUIRED)
    message("Including Boost header directory:\t ${Boost_INCLUDE_DIRS}/boost_1_72_0" )
    include_directories("${Boost_INCLUDE_DIRS}")
endif()

# Everything beyond this point should not need any specifically typed paths...
########################################################################################################################
########################################################################################################################
TARGET_LINK_LIBRARIES(${PROJECT_NAME} PUBLIC ${Boost_LIBRARIES})
set(SOFA_ROOT ${sofa_rooot})
set(SOFA_DIR ${sofa_diiir})

find_package(EIGEN3 REQUIRED)

# Define some macros to make it easier to add each sofa item because it does weird things

macro(add_sofa_core plugin)
    # reset sofa root because find_package keeps unsetting it... so I'm not sure if it finds the right sofa if you have
    # multiple installations...
    set(SOFA_ROOT ${sofa_rooot})
    set(SOFA_DIR ${sofa_diiir})
    find_package(${plugin} REQUIRED)
    set(SOFA_ROOT ${sofa_rooot})
    set(SOFA_DIR ${sofa_diiir})
    set("${plugin}_INCLUDE_DIR" "${SOFA_ROOT}/${plugin}")
    message("including:\t     ${${plugin}_INCLUDE_DIR}")
    include_directories(${${plugin}_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${plugin})
endmacro()

macro(add_sofa_plugin plugin)
    # reset sofa root because find_package keeps unsetting it... so I'm not sure if it finds the right sofa if you have
    # multiple installations...    set(SOFA_ROOT ${sofa_rooot})
    set(SOFA_DIR ${sofa_diiir})
    set("${plugin}_DIR"   "${SOFA_ROOT}/plugins/${plugin}/lib/cmake/${plugin}")
    find_package(${plugin} REQUIRED)
    set(SOFA_ROOT ${sofa_rooot})
    set(SOFA_DIR ${sofa_diiir})
    set("${plugin}_INCLUDE_DIR" "${SOFA_ROOT}/plugins/${plugin}/include")
    message("including:\t     ${${plugin}_INCLUDE_DIR}")
    include_directories(${${plugin}_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${plugin})
endmacro()

# Add all the required sofa stuff
add_sofa_core(SofaFramework)
add_sofa_core(SofaBase)
add_sofa_core(SofaSimulation)
add_sofa_core(SofaCommon)
add_sofa_core(SofaGeneral)
add_sofa_plugin(SofaGraphComponent)
add_sofa_plugin(SofaBoundaryCondition)
add_sofa_plugin(SofaSimpleFem)




