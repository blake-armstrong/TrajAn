add_subdirectory(tacpp/core)
add_subdirectory(tacpp/io)
add_subdirectory(tacpp/main)

set(TACPP_INCLUDE_FILES
    "${TACPP_CORE_INCLUDE_FILES}" "${TACPP_IO_INCLUDE_FILES}"
    "${TACPP_MAIN_INCLUDE_FILES}"
    PARENT_SCOPE)

set(TACPP_INCLUDE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}
    PARENT_SCOPE)
