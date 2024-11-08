cmake_minimum_required(VERSION 3.10)
project(DMPPICFoam VERSION 0.1.0 LANGUAGES C CXX)

# Set the build type
set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_BUILD_TYPE Debug)

# Set optimization levels
if(CMAKE_BUILD_TYPE MATCHES Debug)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")  # No optimization
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
elseif(CMAKE_BUILD_TYPE MATCHES Release)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")  # Full optimization
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
endif()


include(CTest)
enable_testing()

add_executable(DMPPICFoam ./src/DMPPICFoam.C)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

set_target_properties(DMPPICFoam PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

include(CPack)

add_definitions(
-Dlinux64
-DWM_DP
-fpermissive
-DNoRepository
)

include_directories(${PROJECT_SOURCE_DIR}/)
include_directories(${PROJECT_SOURCE_DIR}/include/)
include_directories(${PROJECT_SOURCE_DIR}/include/OSspecific/POSIX/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/lib/DPMTurbulenceModels)
include_directories(${PROJECT_SOURCE_DIR}/include/finiteVolume/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/finiteArea/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/meshTools/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/lagrangian/basic/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/lagrangian/intermediate/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/thermophysicalModels/specie/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/transportModels/compressible/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/thermophysicalModels/basic/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/thermophysicalModels/reactionThermo/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/thermophysicalModels/radiation/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/transportModels)
include_directories(${PROJECT_SOURCE_DIR}/include/transportModels/incompressible/singlePhaseTransportModel)
include_directories(${PROJECT_SOURCE_DIR}/include/TurbulenceModels/turbulenceModels/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/TurbulenceModels/incompressible/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/regionModels/regionModel/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/regionModels/surfaceFilmModels/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/regionFaModels/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/faOptions/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/include/TurbulenceModels/phaseIncompressible/lnInclude)

# include_directories(${PROJECT_SOURCE_DIR}/include/OpenFOAM/db/Time/)
include_directories(${PROJECT_SOURCE_DIR}/include/OpenFOAM/lnInclude/)


target_link_libraries(${PROJECT_NAME}
${PROJECT_SOURCE_DIR}/lib/libfiniteVolume.so
${PROJECT_SOURCE_DIR}/lib/libfvOptions.so
${PROJECT_SOURCE_DIR}/lib/libmeshTools.so
${PROJECT_SOURCE_DIR}/lib/liblagrangian.so
${PROJECT_SOURCE_DIR}/lib/liblagrangianIntermediate.so
${PROJECT_SOURCE_DIR}/lib/liblagrangianTurbulence.so
${PROJECT_SOURCE_DIR}/lib/libspecie.so
${PROJECT_SOURCE_DIR}/lib/libradiationModels.so
${PROJECT_SOURCE_DIR}/lib/libincompressibleTransportModels.so
${PROJECT_SOURCE_DIR}/lib/libturbulenceModels.so
${PROJECT_SOURCE_DIR}/lib/libincompressibleTurbulenceModels.so
${PROJECT_SOURCE_DIR}/lib/libDPMTurbulenceModels.so
${PROJECT_SOURCE_DIR}/lib/libregionModels.so
${PROJECT_SOURCE_DIR}/lib/libsurfaceFilmModels.so
${PROJECT_SOURCE_DIR}/lib/libsampling.so
${PROJECT_SOURCE_DIR}/lib/libregionFaModels.so
${PROJECT_SOURCE_DIR}/lib/libfiniteArea.so
${PROJECT_SOURCE_DIR}/lib/libfaOptions.so
${PROJECT_SOURCE_DIR}/lib/libOpenFOAM.so
${PROJECT_SOURCE_DIR}/lib/sys-openmpi/libPstream.so
${PROJECT_SOURCE_DIR}/lib/libfileFormats.so
${PROJECT_SOURCE_DIR}/lib/libsurfMesh.so
${PROJECT_SOURCE_DIR}/lib/libfluidThermophysicalModels.so
${PROJECT_SOURCE_DIR}/lib/libsolidThermo.so
${PROJECT_SOURCE_DIR}/lib/libcompressibleTurbulenceModels.so
${PROJECT_SOURCE_DIR}/lib/libreactionThermophysicalModels.so
${PROJECT_SOURCE_DIR}/lib/liboverset.so
${PROJECT_SOURCE_DIR}/lib/libdistributionModels.so
${PROJECT_SOURCE_DIR}/lib/libcompressibleTransportModels.so
${PROJECT_SOURCE_DIR}/lib/libthermophysicalProperties.so
${PROJECT_SOURCE_DIR}/lib/libSLGThermo.so
${PROJECT_SOURCE_DIR}/lib/libdynamicMesh.so
${PROJECT_SOURCE_DIR}/lib/libdynamicFvMesh.so
${PROJECT_SOURCE_DIR}/lib/libdistributed.so
${PROJECT_SOURCE_DIR}/lib/libtwoPhaseMixture.so
${PROJECT_SOURCE_DIR}/lib/libthermoTools.so
${PROJECT_SOURCE_DIR}/lib/libsolidSpecie.so
${PROJECT_SOURCE_DIR}/lib/libblockMesh.so
${PROJECT_SOURCE_DIR}/lib/libdecompositionMethods.so
${PROJECT_SOURCE_DIR}/lib/libextrudeModel.so
${PROJECT_SOURCE_DIR}/lib/libincompressibleMultiphaseSystems.so
${PROJECT_SOURCE_DIR}/lib/libgeometricVoF.so
${PROJECT_SOURCE_DIR}/lib/libtwoPhaseProperties.so
${PROJECT_SOURCE_DIR}/lib/libcombustionModels.so
${PROJECT_SOURCE_DIR}/lib/libsaturationModel.so
${PROJECT_SOURCE_DIR}/lib/libinterfaceProperties.so
${PROJECT_SOURCE_DIR}/lib/libchemistryModel.so
${PROJECT_SOURCE_DIR}/lib/libODE.so
${PROJECT_SOURCE_DIR}/lib/libfileFormats.so
)
