

#include <gtest/gtest.h>

#include <parametrics/surfaces/gmpapple.h>
#include <parametrics/surfaces/gmpapple2.h>
#include <parametrics/surfaces/gmpasteroidalsphere.h>
using namespace GMlib;


namespace {

template <typename T>
inline
void testPSurfaceStandardMethodCalls( PSurf<T,3>& surface ) {
  surface.evaluate( surface.getParStartU() + surface.getParDeltaU() * 0.5, 0, surface.getParStartV() + surface.getParDeltaV() * 0.5, 0 );
}


TEST(Parametrics_Surfaces, PAppleCompile) {

    auto psurface = PApple<float>();
    testPSurfaceStandardMethodCalls(psurface);
}

TEST(Parametrics_Surfaces, PApple2Compile) {

    auto psurface = PApple2<float>();
    testPSurfaceStandardMethodCalls(psurface);
}

TEST(Parametrics_Surfaces, PAsteroidalSphereCompile) {

    auto psurface = PAsteroidalSphere<float>();
    testPSurfaceStandardMethodCalls(psurface);
}


}

























