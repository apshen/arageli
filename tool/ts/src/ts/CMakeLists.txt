set(arageli_ts_sources
    ts.cpp
    tstcontainer.cpp
)

add_library(ts STATIC EXCLUDE_FROM_ALL ${arageli_ts_sources})
