﻿cmake_minimum_required (VERSION 3.10)

add_executable (nQuantCpp "nQuantCpp.cpp" "nQuantCpp.h" "nQuantCpp.rc" "bitmapUtilities.cpp" "bitmapUtilities.h" "CIELABConvertor.cpp" "CIELABConvertor.h" "DivQuantizer.cpp" "DivQuantizer.h"
 "Dl3Quantizer.cpp" "Dl3Quantizer.h" "EdgeAwareSQuantizer.cpp" "EdgeAwareSQuantizer.h" "HilbertCurve.cpp" "HilbertCurve.h" "MedianCut.cpp" "MedianCut.h" "MoDEQuantizer.cpp" "MoDEQuantizer.h"
 "NeuQuantizer.cpp" "NeuQuantizer.h" "PnnLABQuantizer.cpp" "PnnLABQuantizer.h" "PnnQuantizer.cpp" "PnnQuantizer.h" "Resource.h"
 "SpatialQuantizer.cpp" "SpatialQuantizer.h" "stdafx.cpp" "stdafx.h"  "WuQuantizer.cpp" "WuQuantizer.h")
target_link_libraries(nQuantCpp gdiplus)
