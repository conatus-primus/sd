//**************************************************************************************//
//                                                                                      //
//    The MIT License (MIT)                                                             //
//                                                                                      //
//    Copyright (c) 2017 ConatusPrimus                                                  //
//                                                                                      //
//    Permission is hereby granted, free of charge, to any person obtaining a copy      //
//    of this software and associated documentation files (the "Software"), to deal     //
//    in the Software without restriction, including without limitation the rights      //
//    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell         //
//    copies of the Software, and to permit persons to whom the Software is             //
//    furnished to do so, subject to the following conditions:                          //
//                                                                                      //
//    The above copyright notice and this permission notice shall be included in all    //
//    copies or substantial portions of the Software.                                   //
//                                                                                      //
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR        //
//    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,          //
//    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
//    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER            //
//    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,     //
//    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE     //
//    SOFTWARE.                                                                         //
//                                                                                      //
//    Sources are located by https://github.com/conatus-primus/SD                       //
//                                                                                      //
//**************************************************************************************//

#include <QDebug>
#include "sdist.h"

#define _USE_MATH_DEFINES 
#include <math.h>

void test1() {
  
  //1. совпадает плоскость х
  double fi[] = {0, M_PI/2, 0};
  double lm[] = {0, 0, M_PI/2 };

  double fi0[] = {0};
  double lm0[] = {0};
  double directRad = {0};

  SD::ambitus a(0);
  for(int i=0; i<sizeof(fi)/sizeof(double); ++i) {    
    a.append(fi[i],lm[i]);
  }
  a.prepare(true);

  double distanceRad = M_PI/180./60.;

  for(int i=0; i<sizeof(fi0)/sizeof(double); ++i) {
    SD::punctumRad g0(fi0[i],lm0[i]);
    SD::punctumRad g1 = SD::operations::toDirection(g0,directRad,distanceRad);

    //перейти в декартовы координаты
    SD::punctumXYZ xyz0 = SD::operations::geo2xyz(g0);
    SD::punctumXYZ xyz1 = SD::operations::geo2xyz(g1);

    //построить плоскость проходяющую через g0 и g1
    //построение плоскости проходящей через центр сферы
    SD::planum inputPlanum = SD::operations::planumCrossZero(xyz0, xyz1);

    std::vector<SD::punctumOnAmbitus> outputPoints;
    a.intersecare(inputPlanum,xyz0,&outputPoints);

    for(unsigned j=0; j<outputPoints.size(); ++j) {
      const SD::punctumOnAmbitus &p = outputPoints[j];
      qDebug() << p.idArcus << p.fi << p.lm;
    }
  }
}

void test2() {
  
  //1. совпадает плоскость х
  double fi[] = {0, M_PI/2, 0};
  double lm[] = {0, 0, M_PI/2 };

  double fi0[] = {0};
  double lm0[] = {M_PI/4};
  double directStep = 5;

  SD::ambitus a(0);
  a.append(fi[1],lm[1]);
  a.append(fi[2],lm[2]);
  a.prepare(false);

  double distanceRad = M_PI/180./60.;

  for(int k=5; k<=360; k+=directStep) {
    double directRad = k/180.*M_PI;

    for(int i=0; i<sizeof(fi0)/sizeof(double); ++i) {

      qDebug() << "--------------------------------";
      qDebug() << "direct=" << k << " {fi,lm}=" <<  fi0[i]/M_PI*180. << "," << lm0[i]/M_PI*180.;

      std::vector<SD::punctumOnAmbitus> outputPoints;
      
      int numberPunctumByDirection;
      a.intersecare(fi0[i],lm0[i],directRad,&outputPoints,numberPunctumByDirection);
      for(unsigned j=0; j<outputPoints.size(); ++j) {
        const SD::punctumOnAmbitus &p = outputPoints[j];
        qDebug() << p.idArcus << ") " << p.fi/M_PI*180 << p.lm/M_PI*180;
      }
    }
  }
}

int main(int argc, char *argv[])
{
  //плоскость проходит через ось у
  SD::punctumRad r0(M_PI_4,M_PI_2);//-M_PI/6.);
  SD::punctumXYZ l0 = SD::operations::geo2xyz(r0);

  SD::punctumRad r1 = SD::operations::toDirection(r0,0,M_PI/180.);
  SD::punctumXYZ l1 = SD::operations::geo2xyz(r1);
  //SD::punctumXYZ l1(1,0,0);

  /*
  SD::planum p01 = SD::operations::planumCrossZero(l0,l1);
  double alfa, beta, gama;
  p01.angles_eyler(alfa, beta, gama);
  qDebug() << SD::operations::G(alfa) << SD::operations::G(beta);
  
  SD::punctumXYZ xyz1 = p01.rotateZ(l0,alfa);
  qDebug() << xyz1.x << xyz1.y << xyz1.z;

  SD::punctumXYZ xyz2 = p01.rotateX(xyz1,beta);
  qDebug() << xyz2.x << xyz2.y << xyz2.z;
  /*
  SD::planum v1;
  v1.a=-1;
  v1.b=0;
  v1.c=0;
  double alfa, beta, gama;
  v1.angles_eyler(alfa, beta, gama);
  qDebug() << SD::operations::G(alfa);
  SD::punctumXYZ xyz; xyz.x=0;xyz.y=1;xyz.z=0;
  SD::punctumXYZ xyz1 = v1.rotate(xyz,alfa,0,0);
  
  {
  SD::planum v1;
  v1.a=1;
  v1.b=1;
  v1.c=0;
  double alfa, beta, gama;
  v1.angles_eyler(alfa, beta, gama);
  qDebug() << SD::operations::G(alfa);
  }

  SD::punctumRad p1(M_PI/4,M_PI/4);
  SD::punctumRad p2(0,0);

  SD::planum pp = SD::operations::planumCrossZero(
    SD::operations::geo2xyz(p1),
    SD::operations::geo2xyz(p2)
  );
  double a1 = SD::operations::anglePlanumsRadWithZ0(pp);


  bool r1 = pp.contains(SD::operations::geo2xyz(p1), SD::operations::EPS);
  r1 = pp.contains(SD::operations::geo2xyz(p2), SD::operations::EPS);

  SD::planumRotate rr(pp); 

  SD::punctumXYZ xyz01 = SD::operations::geo2xyz(p1);
  SD::punctumXYZ xyz02 = SD::operations::geo2xyz(p2);

  SD::punctumXYZ xyz1 = pp.to(xyz01,a1);
  SD::punctumXYZ xyz2 = pp.to(xyz02,a1);
  */
  /*
  SD::punctumXYZ xyz1 = rr.rotateY(xyz01);
  SD::punctumXYZ xyz2 = rr.rotateY(xyz02);

  r1 = rr.contains(SD::operations::geo2xyz(p1), SD::operations::EPS);
  r1 = rr.contains(SD::operations::geo2xyz(p2), SD::operations::EPS);
  */
  /*
  {
  SD::planum p;
  p.a = 1;
  p.b = 1;
  p.c = 0;
  double alfa, beta, gamma;
  p.angles(alfa, beta, gamma);
  qDebug() << SD::operations::G(alfa) << SD::operations::G(beta) << SD::operations::G(gamma);
  
  SD::planumRotate pp;
  pp.a = 1; pp.b = 0; pp.c = 0;
  pp.angles();

  std::vector<SD::punctumXYZ> vxyz;
  for(int i=0; i<360; i+=10) {
    SD::punctumXYZ xyz0, xyz;
    xyz0.x=0;
    xyz0.y = cos(i/180.*M_PI);
    xyz0.z = sin(i/180.*M_PI);

    xyz.z = xyz0.z;
    xyz.x = xyz0.y*cos(M_PI/2.);
    xyz.y = xyz0.y*cos(M_PI/2.);

    if(!pp.contains(xyz, SD::operations::EPS)) {
      continue;
    }
    //SD::punctumXYZ xyz1 = pp.rotateY(xyz);
    SD::punctumXYZ xyz2 = pp.rotateZ(xyz);
    //SD::punctumXYZ xyz3 = pp.rotateX(xyz2);
    qDebug() << xyz2.x << xyz2.y << xyz2.z;
  }
  }*/
  test1();
  test2();
  return 0;
}
