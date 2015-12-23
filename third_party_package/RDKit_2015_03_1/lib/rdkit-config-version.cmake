set (PACKAGE_VERSION "1.2015.03.1")
if (NOT "${PACKAGE_FIND_VERSION}" VERSION_GREATER  "1.2015.03.1")
    set (PACKAGE_VERSION_COMPATIBLE 1) # assuming backward compatible
    if ("$PACKAGE_FIND_VERSION}" VERSION_EQUAL  "1.2015.03.1")
        set (PACKAGE_VERSION_EXACT 1) # matches exactly
    endif ()
endif ()
