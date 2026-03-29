cmake_minimum_required(VERSION 3.10)

project(mmd VERSION 0.1
        DESCRIPTION "MD in modern Fortran"
        LANGUAGES Fortran)
enable_language(Fortran)

set(FVERSION "-std=f2023")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION}")

# Source code
add_subdirectory(src)
