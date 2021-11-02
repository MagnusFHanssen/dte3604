
#include <iostream>

#include "scenario.h"
#include "testtorus.h"
#include "custom/lotus.h"
#include "custom/bsplinecustom.h"


// hidmanager
#include "hidmanager/defaulthidmanager.h"

// gmlib
#include <scene/light/gmpointlight.h>
#include <scene/sceneobjects/gmpathtrack.h>
#include <scene/sceneobjects/gmpathtrackarrows.h>

// qt
#include <QQuickItem>


template <typename T>
inline
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
  out << v.size() << std::endl;
  for(unsigned int i=0; i<v.size(); i++) out << " " << v[i];
  out << std::endl;
  return out;
}




void Scenario::initializeScenario() {

  // Insert a light
  GMlib::Point<GLfloat,3> init_light_pos( 2.0, 4.0, 10 );
  GMlib::PointLight *light = new GMlib::PointLight(  GMlib::GMcolor::white(), GMlib::GMcolor::white(),
                                                     GMlib::GMcolor::white(), init_light_pos );
  light->setAttenuation(0.8f, 0.002f, 0.0008f);
  this->scene()->insertLight( light, false );

  // Insert Sun
  this->scene()->insertSun();

  // Default camera parameters
  int init_viewport_size = 600;
  GMlib::Point<float,3>  init_cam_pos( 0.0f, 0.0f, 0.0f );
  GMlib::Vector<float,3> init_cam_dir( 0.0f, 1.0f, 0.0f );
  GMlib::Vector<float,3> init_cam_up(  1.0f, 0.0f, 0.0f );

  // Projection cam
  auto proj_rcpair = createRCPair("Projection");
  proj_rcpair.camera->set(init_cam_pos,init_cam_dir,init_cam_up);
  proj_rcpair.camera->setCuttingPlanes( 1.0f, 8000.0f );
  proj_rcpair.camera->rotateGlobal( GMlib::Angle(-45), GMlib::Vector<float,3>( 1.0f, 0.0f, 0.0f ) );
  proj_rcpair.camera->translateGlobal( GMlib::Vector<float,3>( 0.0f, -20.0f, 20.0f ) );
  scene()->insertCamera( proj_rcpair.camera.get() );
  proj_rcpair.renderer->reshape( GMlib::Vector<int,2>(init_viewport_size, init_viewport_size) );


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

//  auto lotus = new Custom::Lotus();
//  lotus->toggleDefaultVisualizer();
//  lotus->sample(400);
//  this->scene()->insert(lotus);
//  lotus->setColor(GMlib::Color(153, 255 ,255));

  GMlib::DVector<GMlib::Point<double,3>> controlPoints(9);
  controlPoints[0] = GMlib::Point<double,3>(0, 0, 0);
  controlPoints[1] = GMlib::Point<double,3>(0, 1, 0);
  controlPoints[2] = GMlib::Point<double,3>(3, 0, 1);
  controlPoints[3] = GMlib::Point<double,3>(3, 4, 0);
  controlPoints[4] = GMlib::Point<double,3>(0, 1, 6);
  controlPoints[5] = GMlib::Point<double,3>(1, 1, 1);
  controlPoints[6] = GMlib::Point<double,3>(0, 8, 0);
  controlPoints[7] = GMlib::Point<double,3>(0, 1, 8);
  controlPoints[8] = GMlib::Point<double,3>(9, 0, 1);


  auto spline = new Custom::BSplineCustom<double>(controlPoints);
  spline->toggleDefaultVisualizer();
  spline->sample(200, 0);
  this->scene()->insert(spline);
  spline->setColor(GMlib::GMcolor::red());

  auto spline2 = new Custom::BSplineCustom<double>(controlPoints);
  spline2->toggleDefaultVisualizer();
  spline2->setBlending(true);
  spline2->sample(200, 0);
  this->scene()->insert(spline2);
  spline2->setColor(GMlib::GMcolor::green());
}




void Scenario::cleanupScenario() {

}




void Scenario::callDefferedGL() {

  GMlib::Array< const GMlib::SceneObject*> e_obj;
  this->scene()->getEditedObjects(e_obj);

  for(int i=0; i < e_obj.getSize(); i++)
    if(e_obj(i)->isVisible()) e_obj[i]->replot();
}

