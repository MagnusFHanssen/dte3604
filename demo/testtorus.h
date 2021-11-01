#ifndef TESTTORUS_H
#define TESTTORUS_H


// gmlib
#include <parametrics/surfaces/gmptorus.h>


class TestTorus : public GMlib::PTorus<float> {
public:
  using PTorus::PTorus;

  ~TestTorus() override {

    if(m_test01)
      remove(test_01_torus.get());
  }

  void test01() {

    GMlib::Vector<float,3> d = evaluate(0.0f,0.0f,0,0)[0][0];
    test_01_torus = std::make_shared<TestTorus,float,float,float>(1.5f,0.5f,0.5f);

    test_01_torus->translate(d + d.getNormalized()*2.0f);
    test_01_torus->rotate( GMlib::Angle(90), GMlib::Vector<float,3>( 0.0f, 1.0f, 0.0f) );
    test_01_torus->toggleDefaultVisualizer();
    test_01_torus->sample(200,200,1,1);
    insert(test_01_torus.get());

    m_test01 = true;
  }


protected:
  void localSimulate(double dt) override {

      static double t=0;
      t += dt;

      GMlib::Vector<float,3> vec(sin(t),cos(t),5*dt);
      vec *= 0.05;
      this->move(vec);
  }

private:
  bool m_test01 {false};
  std::shared_ptr<TestTorus> test_01_torus {nullptr};

}; // END class TestTorus



#endif // TESTTORUS_H
