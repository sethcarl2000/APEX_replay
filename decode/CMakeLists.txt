
#-------------------------------------------------------------------------------------------
#
#   First, we're going to add the 'decode' libaray, which handles operations with CODA files
#
set(src_decode
    count_physics_events.cxx
    init_decoders.cxx
    decode_coda_file.cxx 
#    Shower.cxx
)

set(headers_decode
    include/APEX/decode.h
#    include/APEX/decode/Shower.h
)

add_library(APEX_decode SHARED ${src_decode})

target_include_directories(APEX_decode PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

# link it with Podd
target_link_libraries(APEX_decode PUBLIC Podd::Podd Podd::HallA APEX_utils)

