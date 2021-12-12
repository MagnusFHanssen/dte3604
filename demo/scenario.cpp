
#include <iostream>

#include "scenario.h"
#include "testtorus.h"
#include "custom/lotus.h"
#include "custom/bsplinecustom.h"
#include "custom/laneriesenfeldclosed.h"
#include "custom/blendingspline.h"
#include "custom/modelsurface.h"
#include "custom/ctorus.h"
#include "custom/blendingsurface.h"


// hidmanager
#include "hidmanager/defaulthidmanager.h"

// gmlib
#include <scene/light/gmpointlight.h>
#include <scene/sceneobjects/gmpointlightg.h>
#include <scene/sceneobjects/gmspotlightg.h>
#include <scene/sceneobjects/gmpathtrack.h>
#include <scene/sceneobjects/gmpathtrackarrows.h>

#include <parametrics/visualizers/gmpsurfnormalsvisualizer.h>
#include <parametrics/visualizers/gmpsurfderivativesvisualizer.h>

#include <parametrics/surfaces/gmpsphere.h>
#include <parametrics/curves/gmpline.h>
#include <parametrics/visualizers/gmpcurvepointsvisualizer.h>

#include <parametrics/curves/gmpcircle.h>

// qt
#include <QQuickItem>

template <typename T>
inline std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  out << v.size() << std::endl;
  for (unsigned int i = 0; i < v.size(); i++)
    out << " " << v[i];
  out << std::endl;
  return out;
}

void Scenario::initializeScenario() {

  // Insert a light
  GMlib::Point<GLfloat, 3> init_light_pos(2.0, 4.0, 10);
  GMlib::PointLightG *light =
      new GMlib::PointLightG(GMlib::GMcolor::white(), GMlib::GMcolor::white(),
                            GMlib::GMcolor::white(), init_light_pos);
  light->setAttenuation(0.8f, 0.002f, 0.0008f);
  this->scene()->insertLight(light, false);
  //this->scene()->insert(light);

  // Insert a spotlight
  GMlib::Point<float, 3> init_spotlight_pos(2.0, 4.0, 10);
  GMlib::Vector<float, 3> spotlight_direction(0, -1.0, 0);
  GMlib::SpotLightG *spotlight =
      new GMlib::SpotLightG(GMlib::GMcolor::white(), GMlib::GMcolor::white(),
                            GMlib::GMcolor::white(), init_spotlight_pos, spotlight_direction, GMlib::Angle(20), 2.0);
  spotlight->setAttenuation(0.8f, 0.002f, 0.0008f);
  this->scene()->insertLight(spotlight, false);
  //this->scene()->insert(spotlight);

  // Insert Sun
  this->scene()->insertSun();

  // Default camera parameters
  int init_viewport_size = 600;
  GMlib::Point<float, 3> init_cam_pos(0.0f, 0.0f, 0.0f);
  GMlib::Vector<float, 3> init_cam_dir(0.0f, 1.0f, 0.0f);
  GMlib::Vector<float, 3> init_cam_up(1.0f, 0.0f, 0.0f);

  // Projection cam
  auto proj_rcpair = createRCPair("Projection");
  proj_rcpair.camera->set(init_cam_pos, init_cam_dir, init_cam_up);
  proj_rcpair.camera->setCuttingPlanes(1.0f, 8000.0f);
  proj_rcpair.camera->rotateGlobal(GMlib::Angle(-45),
                                   GMlib::Vector<float, 3>(1.0f, 0.0f, 0.0f));
  proj_rcpair.camera->translateGlobal(
      GMlib::Vector<float, 3>(0.0f, -20.0f, 20.0f));
  scene()->insertCamera(proj_rcpair.camera.get());
  proj_rcpair.renderer->reshape(
      GMlib::Vector<int, 2>(init_viewport_size, init_viewport_size));

  //proj_rcpair.renderer->setClearColor(GMlib::GMcolor::lightCyan());

  /***************************************************************************
   *                                                                         *
   * Standar example, including path track and path track arrows             *
   *                                                                         *
   ***************************************************************************/

  //  GMlib::Material mm(GMlib::GMmaterial::polishedBronze());
  //  mm.set(45.0);

  //  auto ptom = new TestTorus(1.0f, 0.4f, 0.6f);
  //  ptom->toggleDefaultVisualizer();
  //  ptom->sample(60,60,1,1);
  //  this->scene()->insert(ptom);
  //  auto ptrack = new GMlib::PathTrack();
  //  ptrack->setLineWidth(2);
  //  ptom->insert(ptrack);
  //  auto ptrack2 = new GMlib::PathTrackArrows();
  //  ptrack2->setArrowLength(2);
  //  ptom->insert(ptrack2);

    auto lotus = new Custom::Lotus();
    lotus->toggleDefaultVisualizer();
    lotus->sample(400);
    //this->scene()->insert(lotus);
    lotus->setColor(GMlib::Color(153, 255 ,255));


  // The controls for quicly swapping between several configuration without comments
  bool showDottedLines = false;
  bool showControlPoly = true;


  GMlib::DVector<GMlib::Point<double, 3>> controlPoints(9);
  controlPoints[0] = GMlib::Point<double, 3>(0, 0, 0);
  controlPoints[1] = GMlib::Point<double, 3>(0, 1, 0);
  controlPoints[2] = GMlib::Point<double, 3>(3, 0, 1);
  controlPoints[3] = GMlib::Point<double, 3>(3, 4, 0);
  controlPoints[4] = GMlib::Point<double, 3>(0, 1, 6);
  controlPoints[5] = GMlib::Point<double, 3>(1, 1, 1);
  controlPoints[6] = GMlib::Point<double, 3>(0, 8, 0);
  controlPoints[7] = GMlib::Point<double, 3>(0, 1, 8);
  controlPoints[8] = GMlib::Point<double, 3>(9, 0, 1);

//  for (int i = 0; i < 9 && showControlPoly; i++) {
//    auto sph = new GMlib::PSphere<double>(0.1);
//    sph->toggleDefaultVisualizer();
//    sph->translate(controlPoints[i]);
//    sph->sample(8, 8, 1, 1);
//    this->scene()->insert(sph);
//  }
//  for (int i = 0; i < 8 && showControlPoly; i++) {
//    auto lin = new GMlib::PLine<double>(controlPoints[i], controlPoints[i + 1]);
//    lin->toggleDefaultVisualizer();
//    lin->sample(2);
//    lin->setColor(GMlib::GMcolor::aliceBlue());
//    this->scene()->insert(lin);
//  }

  auto c_viz = new GMlib::PCurvePointsVisualizer<double, 3>(3);
  auto c_viz2 = new GMlib::PCurvePointsVisualizer<double, 3>(3);
  c_viz2->setColor(GMlib::GMcolor::green());

//  auto spline = new Custom::BSplineCustom<double>(controlPoints);
//  spline->setColor(GMlib::GMcolor::red());
//  if (showDottedLines) {
//    spline->insertVisualizer(c_viz);
//  } else {
//    spline->toggleDefaultVisualizer();
//  }
//  spline->sample(200, 0);
//  this->scene()->insert(spline);


//  auto spline2 = new Custom::BSplineCustom<double>(controlPoints);
//  spline2->setColor(GMlib::GMcolor::green());
//  if (showDottedLines) {
//      spline2->insertVisualizer(c_viz2);
//  } else {
//      spline2->toggleDefaultVisualizer();
//  }
//  spline2->setBlending(true);
//  spline2->sample(200, 0);
//    this->scene()->insert(spline2);

//  GMlib::DVector<GMlib::Point<double,3>> points(12);
//  for (int i  =0; i < 12; i++){
//      points[i] = GMlib::Point<double,3>(10*cos(i * M_PI/6),
//                                         10*sin(i * M_PI/6),
//                                          0.0);
//  }

//  auto spline3 = new Custom::BSplineCustom<double>(points, 5);
//  spline3->setColor(GMlib::GMcolor::yellow());
//  spline3->toggleDefaultVisualizer();
//  spline3->setBlending(false);
//  spline3->sample(200, 0);
//  this->scene()->insert(spline3);

//  auto lr = new Custom::LaneRiesenfeldClosed<double>(controlPoints, 2);
//  lr->sample(4, 0);
//  lr->toggleDefaultVisualizer();
//  this->scene()->insert(lr);

//  auto circle = new GMlib::PCircle<double>();

  auto bls = new Custom::BlendingSpline<double>(lotus, 10);
  bls->toggleDefaultVisualizer();
  //bls->insertVisualizer(c_viz);
  bls->sample(200, 0);
  //bls->showControlCurves();
  //this->scene()->insert(bls);

  auto surf1 = new Custom::CTorus<float>(); //new Custom::ModelSurface<double>(3, 1.0, 0);

  //auto n_vis = new GMlib::PSurfNormalsVisualizer<double, 3>();
  //auto d_vis = new GMlib::PSurfDerivativesVisualizer<double, 3>(0, 1);

  surf1->toggleDefaultVisualizer();
  //surf1->insertVisualizer(d_vis);

  surf1->sample(20, 20, 1, 1);
  this->scene()->insert(surf1);
  surf1->setMaterial(GMlib::GMmaterial::emerald());

  auto b_surf = new Custom::BlendingSurface<float>(surf1, 10, 10);
  b_surf->toggleDefaultVisualizer();
  b_surf->sample(20, 20);
  b_surf->showLocalSurfaces();
  this->scene()->insert(b_surf);
}

void Scenario::cleanupScenario() {}

void Scenario::callDefferedGL() {

  GMlib::Array<const GMlib::SceneObject *> e_obj;
  this->scene()->getEditedObjects(e_obj);

  for (int i = 0; i < e_obj.getSize(); i++)
    if (e_obj(i)->isVisible())
      e_obj[i]->replot();
}
