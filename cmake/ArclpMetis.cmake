# Optional nested-dissection ordering through METIS.  When found, the
# analysis also tries it and keeps whichever ordering gives the smaller
# factor.  Sets ARCLP_USE_METIS, ARCLP_METIS_INCLUDE_DIR, ARCLP_METIS_LIBRARY.

set(ARCLP_USE_METIS OFF)
if(NOT ARCLP_METIS STREQUAL "off")
    find_path(ARCLP_METIS_INCLUDE_DIR metis.h)
    find_library(ARCLP_METIS_LIBRARY metis)
    if(ARCLP_METIS_INCLUDE_DIR AND ARCLP_METIS_LIBRARY)
        set(ARCLP_USE_METIS ON)
    elseif(ARCLP_METIS STREQUAL "on")
        message(FATAL_ERROR "ARCLP_METIS=on but METIS was not found")
    endif()
endif()
message(STATUS "arclp METIS ordering: ${ARCLP_USE_METIS}")
