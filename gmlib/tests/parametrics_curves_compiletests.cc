

#include <gtest/gtest.h>

#include <parametrics/curves/gmpcircle.h>
#include <parametrics/curves/gmpbutterfly.h>
#include <parametrics/curves/gmprosecurve.h>
#include <parametrics/curves/gmpchrysanthemumcurve.h>
#include <parametrics/curves/gmpbasiscurve.h>
#include <parametrics/curves/gmparc.h>
#include <parametrics/curves/gmpbeziercurve.h>
#include <parametrics/curves/gmpbsplinebasiscurve.h>
#include <parametrics/curves/gmpbsplinecurve.h>
#include <parametrics/curves/gmpline.h>
#include <parametrics/curves/gmpsubcurve.h>
#include <parametrics/surfaces/gmpapple.h>
#include <parametrics/curves/gmpsurfcurve.h>
#include <parametrics/triangles/gmpbeziertriangle.h>
#include <parametrics/curves/gmptriangcurve.h>
using namespace GMlib;


namespace {

template <typename T>
inline
void testPCurveStandardMethodCalls( PCurve<T,3>& curve ) {
  curve.evaluate( curve.getParStart() + curve.getParDelta() * 0.5, 0 );
}


TEST(Parametrics_Curves, PCircleCompile) {

  auto pcircle = PCircle<float>();
  testPCurveStandardMethodCalls(pcircle);
}

TEST(Parametrics_Curves, PButterflyCompile) {

  auto pbutterfly = PButterfly<float>();
  testPCurveStandardMethodCalls(pbutterfly);
}

TEST(Parametrics_Curves, PRoseCurveCompile) {

  auto prose = PRoseCurve<float>();
  testPCurveStandardMethodCalls(prose);
}

TEST(Parametrics_Curves, PChrysanthemumCurveCompile) {

  auto pchrysanthemum = PChrysanthemumCurve<float>();
  testPCurveStandardMethodCalls(pchrysanthemum);
}

TEST(Parametrics_Curves, PBasisCurveCompile) {

  auto pbasis = PBasisCurve<float>();
  testPCurveStandardMethodCalls(pbasis);
}

TEST(Parametrics_Curves, PArcCompile) {

  float speed = 1;
  float curvature = 1;
  float start = 0;
  float end = 3.14;
  GMlib::DVector<GMlib::Vector<float,3>> vec;
  vec.setDim(3);
  vec[0]=GMlib::Vector<float,3>(0,0,0);
  vec[1]=GMlib::Vector<float,3>(1,1,0);
  vec[2]=GMlib::Vector<float,3>(2,0,0);

  auto parc1 = PArc<float>(speed,curvature);
  auto parc2 = PArc<float>(speed, curvature, start, end);
  auto parc3 = PArc<float>(vec, start, 0.5, end);
  testPCurveStandardMethodCalls(parc1);
  testPCurveStandardMethodCalls(parc2);
  testPCurveStandardMethodCalls(parc3);
}

TEST(Parametrics_Curves, PBezierCurveCompile) {

  GMlib::DVector<GMlib::Vector<float,3>> vec;
  vec.setDim(3);
  vec[0]=GMlib::Vector<float,3>(0,0,0);
  vec[1]=GMlib::Vector<float,3>(1,1,0);
  vec[2]=GMlib::Vector<float,3>(2,0,0);

  auto pbezier1 = PBezierCurve<float>(vec);
  auto pbezier2 = PBezierCurve<float>(vec, 0.0, 0.5, 1.0);
  testPCurveStandardMethodCalls(pbezier1);
  testPCurveStandardMethodCalls(pbezier2);
}

TEST(Parametrics_Curves, PBSplineBasisCurveCompile) {

  GMlib::DVector<float> vec;
  vec.setDim(5);
  vec[0]=1.0;
  vec[1]=2.0;
  vec[2]=3.0;
  vec[3]=4.0;
  vec[4]=5.0;

  auto psplinebasis = PBSplineBasisCurve<float>(vec);
  testPCurveStandardMethodCalls(psplinebasis);
}

TEST(Parametrics_Curves, PBSplineCurveCompile) {

  GMlib::DVector<GMlib::Vector<float,3>> vec;
  vec.setDim(6);
  vec[0]=GMlib::Vector<float,3>(0,0,0);
  vec[1]=GMlib::Vector<float,3>(1,1,0);
  vec[2]=GMlib::Vector<float,3>(2,0,0);
  vec[3]=GMlib::Vector<float,3>(2,1,0);
  vec[4]=GMlib::Vector<float,3>(3,1,0);
  vec[5]=GMlib::Vector<float,3>(4,0,0);

  auto pbspline = PBSplineCurve<float>(vec, 1);
  testPCurveStandardMethodCalls(pbspline);
}

TEST(Parametrics_Curves, PLineCompile) {

  GMlib::Point<float,3> p1 = GMlib::Point<float,3>(0,0,0);
  GMlib::Point<float,3> p2 = GMlib::Point<float,3>(0,0,2);

  auto pline1 = PLine<float>(p1, GMlib::Vector<float,3>(1,1,0));
  auto pline2 = PLine<float>(p1, p2);
  testPCurveStandardMethodCalls(pline1);
  testPCurveStandardMethodCalls(pline2);
}

TEST(Parametrics_Curves, PSubCurveCompile) {

  auto pcircle = new PCircle<float>();
  float start = 0.0;
  float end = 3.14;

  auto psub1 = PSubCurve<float>(pcircle, start, end);
  auto psub2 = PSubCurve<float>(pcircle, start, end, 1.57f);
  testPCurveStandardMethodCalls(psub1);
  testPCurveStandardMethodCalls(psub2);
}

TEST(Parametrics_Curves, PSurfCurveCompile) {

  auto papple = new PApple<float>();
  GMlib::Point<float,2> p1 = GMlib::Point<float,2>(0,0);
  GMlib::Point<float,2> p2 = GMlib::Point<float,2>(0,2);

  auto psurf1 = PSurfCurve<float>(papple, p1, p2);
  auto psurf2 = PSurfCurve<float>(papple, p1, p2, GMlib::Vector<float,2>(1,0), GMlib::Vector<float,2>(-1,0));
  testPCurveStandardMethodCalls(psurf1);
  testPCurveStandardMethodCalls(psurf2);
}

//TEST(Parametrics_Curves, PTriangCurveCompile) {

//  GMlib::DVector<GMlib::Vector<float,3>> vec;
//  vec.setDim(3);
//  vec[0]=GMlib::Vector<float,3>(0,0,0);
//  vec[1]=GMlib::Vector<float,3>(1,1,0);
//  vec[2]=GMlib::Vector<float,3>(2,0,0);
//  auto pbeziertri = new PBezierTriangle<float>(vec);
//  GMlib::Point<float,3> p1 = GMlib::Point<float,3>(1,0,0);
//  GMlib::Point<float,3> p2 = GMlib::Point<float,3>(0,0,1);

//  auto ptri1 = PTriangCurve<float>(pbeziertri, p1, p2);
//  auto ptri2 = PTriangCurve<float>(pbeziertri, p1, p2, GMlib::Vector<float,3>(0,1,0), GMlib::Vector<float,3>(0,-1,0));
//  testPCurveStandardMethodCalls(ptri1);
//  testPCurveStandardMethodCalls(ptri2);
//}

}

























