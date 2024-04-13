include_directories         (${CMAKE_SOURCE_DIR}/src)

# add_executable              (vtk_volume                       vtk_volume.cpp)

# target_link_libraries       (vtk_volume                       ${libraries})

add_executable(fixed fixed.cpp)
add_executable(multiblock_noblend_continuous multiblock_noblend_continuous.cpp)
target_link_libraries(fixed  ${libraries})
target_link_libraries(multiblock_noblend_continuous ${libraries})

set_property                (TARGET fixed                PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
set_property                (TARGET multiblock_noblend_continuous  PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

install                     (TARGETS fixed multiblock_noblend_continuous
                            DESTINATION ${CMAKE_INSTALL_PREFIX}/src/fixed/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)





