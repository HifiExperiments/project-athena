macro(TARGET_ASSIMP)

	if(WIN32)
	    set(VCVER vc140 vc141 vc142 )
	    set(CRT mt md)
	    set(DBG_NAMES)
	    set(REL_NAMES)
	    foreach(_ver IN LISTS VCVER)
	        foreach(_crt IN LISTS CRT)
	            list(APPEND DBG_NAMES assimp-${_ver}-${_crt}d)
	            list(APPEND REL_NAMES assimp-${_ver}-${_crt})
	        endforeach()
	    endforeach()
	endif()

    find_library(ASSIMP_LIBRARY_RELEASE assimp ${REL_NAMES} PATHS ${VCPKG_INSTALL_ROOT}/lib NO_DEFAULT_PATH)
    find_library(ASSIMP_LIBRARY_RELEASE assimp ${DBG_NAMES} PATHS ${VCPKG_INSTALL_ROOT}/lib NO_DEFAULT_PATH)
    select_library_configurations(ASSIMP)

    target_link_libraries(${TARGET_NAME} ${ASSIMP_LIBRARIES})

endmacro()
