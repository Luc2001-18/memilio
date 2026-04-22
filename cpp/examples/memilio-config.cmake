get_filename_component(MEMILIO_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

if(NOT TARGET memilio::memilio)
    include("${MEMILIO_CMAKE_DIR}/memilio-targets.cmake")
endif()

set(MEMILIO_LIBRARIES memilio::memilio)
