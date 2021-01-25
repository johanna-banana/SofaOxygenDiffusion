cmake_minimum_required(VERSION 3.17)
project(SofaOxygenDiffusion)
set(CPP_FILE_NAMES
        SofaOxygenDiffusion.cpp SofaOxygenDiffusion.h
        TetrahedronO2DiffusionFEMForceField.inl TetrahedronO2DiffusionFEMForceField.h)
set(test_executable ${PROJECT_NAME}_test)
add_executable(${test_executable} test.cpp ${CPP_FILE_NAMES})
target_compile_features(${test_executable} PRIVATE cxx_std_17)

set(MAKE_MEX_FILE false)  # set to true to compile Matlab interface
set(compile_with_gui false)  # Not working yet. leave false.

include_directories(./)
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/install)

# Define needed core items
set(sofa_core_items)
list(APPEND sofa_core_items SofaFramework)
list(APPEND sofa_core_items SofaBase)
list(APPEND sofa_core_items SofaSimulation)
list(APPEND sofa_core_items SofaCommon)
list(APPEND sofa_core_items SofaGeneral)

# Define needed plugins
set(sofa_plugin_items)
list(APPEND sofa_plugin_items SofaGraphComponent)
list(APPEND sofa_plugin_items SofaBoundaryCondition)
list(APPEND sofa_plugin_items SofaSimpleFem)
list(APPEND sofa_plugin_items SofaGeneralLoader)
list(APPEND sofa_plugin_items SofaImplicitOdeSolver)

if(compile_with_gui)
    list(APPEND sofa_core_items SofaGui)
    list(APPEND sofa_plugin_items SofaOpenglVisual)
endif()

if (MAKE_MEX_FILE)
    include(FindMatlab)
    SET( CMAKE_CXX_COMPILER mex )
    SET( CMAKE_C_COMPILER   mex )
    set(MATLAB_TARGET ${PROJECT_NAME})
    matlab_add_mex(NAME ${MATLAB_TARGET} SRC ${CPP_FILE_NAMES})
    target_compile_features(${MATLAB_TARGET} PRIVATE cxx_std_17)
    set(TARGET_SUFFIX mexw64)
    set_target_properties(${MATLAB_TARGET} PROPERTIES SUFFIX ".${TARGET_SUFFIX}")
    message(STATUS ${Matlab_ROOT_DIR})
    INCLUDE_DIRECTORIES("${Matlab_ROOT_DIR}/extern/include")
    target_compile_definitions(${MATLAB_TARGET} PRIVATE MAKE_MEX_FILE)
else()
    set(CMAKE_CXX_STANDARD ++17)
endif ()

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
#TARGET_LINK_LIBRARIES(${PROJECT_NAME} PUBLIC ${Boost_LIBRARIES})
set(SOFA_ROOT ${sofa_rooot})
set(SOFA_DIR ${sofa_diiir})

if(compile_with_gui)
    add_compile_definitions(USE_GUI)
endif()

find_package(EIGEN3 REQUIRED)

# Define some macros to make it easier to add each sofa item because it does weird things
set(MY_DEPENDENCY_PATHS "")
macro(add_sofa_core target plugin)
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
    target_link_libraries(${target} ${plugin})
    list(APPEND MY_DEPENDENCY_PATHS "${SOFA_ROOT}/${plugin}/bin")
endmacro()

macro(add_sofa_plugin target plugin)
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
    target_link_libraries(${target} ${plugin})
    list(APPEND MY_DEPENDENCY_PATHS "${SOFA_ROOT}/plugins/${plugin}/bin")
endmacro()


# Add all the required sofa stuff
foreach(_core_item ${sofa_core_items})
    add_sofa_core(${test_executable} ${_core_item})
endforeach()
foreach(_plugin_item ${sofa_plugin_items})
    add_sofa_plugin(${test_executable} ${_plugin_item})
endforeach()

if (MAKE_MEX_FILE)
    find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY ENG_LIBRARY MEX_COMPILER )
    target_link_libraries(${MATLAB_TARGET} ${MATLAB_LIBRARIES})

    target_compile_definitions(${MATLAB_TARGET} PUBLIC
            MATLAB_MEX_FILE
            _CRT_SECURE_NO_WARNINGS)

    # Add all the required sofa stuff
    foreach(_core_item ${sofa_core_items})
        add_sofa_core(${MATLAB_TARGET} ${_core_item})
    endforeach()
    foreach(_plugin_item ${sofa_plugin_items})
        add_sofa_plugin(${MATLAB_TARGET} ${_plugin_item})
    endforeach()

    # Transfer the value of ${MY_DEPENDENCY_PATHS} into the install script
    install(CODE "set(MY_DEPENDENCY_PATHS \"${MY_DEPENDENCY_PATHS}\")")
    install(CODE "set(mylibtarget \"${MATLAB_TARGET}\")")
    install(CODE [[
      set(CONFIGURATIONS Release)
      file(GET_RUNTIME_DEPENDENCIES
        LIBRARIES   #$<TARGET_FILE:SofaOxygenDiffusion>
        RESOLVED_DEPENDENCIES_VAR _r_deps
        UNRESOLVED_DEPENDENCIES_VAR _u_deps
        DIRECTORIES ${MY_DEPENDENCY_PATHS}
      )
      foreach(_file ${_r_deps})
        file(INSTALL
          DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
          TYPE SHARED_LIBRARY
          FOLLOW_SYMLINK_CHAIN
          FILES "${_file}"
        )
      endforeach()
      list(LENGTH _u_deps _u_length)
      if("${_u_length}" GREATER 0)
        message(WARNING "${_u_length} Unresolved dependencies detected!")
      endif()
    ]])
endif ()