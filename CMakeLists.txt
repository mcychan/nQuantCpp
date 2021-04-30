# CMakeList.txt: 最上層的 CMake 專案檔，執行全域設定
# 並在這裡包含子專案。
#
cmake_minimum_required (VERSION 3.10)

project ("nQuantCpp")
add_definitions(-DUNICODE -D_UNICODE)
find_library(GDIPLUS_LIBRARY NAMES libgdiplus gdiplus)
set(CMAKE_CXX_STANDARD 17)
set(GDIPLUS_LIBRARY gdiplus)

# 包含子專案。
add_subdirectory ("nQuantCpp")
