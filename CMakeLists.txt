cmake_minimum_required (VERSION 3.10)

project ("nQuantCpp")
add_definitions(-DUNICODE -D_UNICODE)
set(CMAKE_CXX_STANDARD 17)
find_file(GDIPLUS_LIBRARY NAMES libgdiplus gdiplus)
if(NOT EXISTS "${GDIPLUS_LIBRARY}")
	find_library(GDIPLUS_LIBRARY NAMES libgdiplus gdiplus)
endif()
message(gdiplus="${GDIPLUS_LIBRARY}")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE
        "MinSizeRel"
        CACHE STRING "Build type: Debug, Release, RelWithDebInfo or MinSizeRel"
              FORCE)

add_subdirectory ("nQuantCpp")
