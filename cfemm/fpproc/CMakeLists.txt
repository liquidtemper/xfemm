add_library(fpproc STATIC
    fpproc.cpp
    makemask.cpp
    CMPointVals.cpp
    CPostProcMElement.cpp
    )
target_include_directories(fpproc PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
target_link_libraries(fpproc PUBLIC femm)

add_executable(fpproc-test
    main.cpp
    )

target_link_libraries(fpproc-test fpproc)

add_executable(fpproc-vgap
    vgap_benchmark.cpp
    )

target_link_libraries(fpproc-vgap fpproc)

install(
    TARGETS fpproc-vgap
    RUNTIME DESTINATION bin
    COMPONENT "cli")
install(
    TARGETS fpproc-test
    RUNTIME DESTINATION bin
    COMPONENT "cli (debug)"
    CONFIGURATIONS "Debug")
# vi:expandtab:tabstop=4 shiftwidth=4:
