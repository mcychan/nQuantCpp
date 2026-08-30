cmake_minimum_required (VERSION 3.10)

if(NOT WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I /usr/share/mingw-w64/include")
endif()

if(MINGW)
  add_compile_options(-municode)
  add_link_options(-municode -Wl,--subsystem,console)
endif()

add_executable(nQuantCpp "nQuantCpp.cpp" "nQuantCpp.h" "nQuantCpp.rc" "bitmapUtilities.cpp" "bitmapUtilities.h" "BlueNoise.cpp" "BlueNoise.h" "CIELABConvertor.cpp" "CIELABConvertor.h" "DivQuantizer.cpp" "DivQuantizer.h"
    "Dl3Quantizer.cpp" "Dl3Quantizer.h" "EdgeAwareSQuantizer.cpp" "EdgeAwareSQuantizer.h" "GifWriter.cpp" "GifWriter.h" "GilbertCurve.cpp" "GilbertCurve.h" "MedianCut.cpp" "MedianCut.h" "Otsu.cpp" "Otsu.h"
    "DblGNGQuantizer.cpp" "DblGNGQuantizer.h" "DblGNGGAQuantizer.cpp" "DblGNGGAQuantizer.h" "NeuQuantizer.cpp" "NeuQuantizer.h" "PnnLABQuantizer.cpp" "PnnLABQuantizer.h" "PnnLABGAQuantizer.cpp" "PnnLABGAQuantizer.h" "PnnQuantizer.cpp" "PnnQuantizer.h" "Resource.h"
    "SpatialQuantizer.cpp" "SpatialQuantizer.h" "stdafx.cpp" "stdafx.h" "WuQuantizer.cpp" "WuQuantizer.h"
    "NsgaIII.cpp" "NsgaIII.h" "APNsgaIII.cpp" "APNsgaIII.h")

target_link_libraries(nQuantCpp PRIVATE gdiplus)
