set(MY_MODS
  nosehoover.f90
  CACHE INTERNAL ""
)
add_library(mylib "${MY_MODS}")
add_executable(mmd main.f90)
target_link_libraries(mmd) PRIVATE mylib)
