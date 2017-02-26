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

  SD::ambitus a;
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

    std::vector<SD::pointInArc> outputPoints;
    a.intersecare(inputPlanum,xyz0,&outputPoints);

    for(unsigned j=0; j<outputPoints.size(); ++j) {
      const SD::pointInArc &p = outputPoints[j];
      qDebug() << p.arcNum << p.geo.fi << p.geo.lm;
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

  SD::ambitus a;
  a.append(fi[1],lm[1]);
  a.append(fi[2],lm[2]);
  a.prepare(false);

  double distanceRad = M_PI/180./60.;

  for(int k=5; k<=360; k+=directStep) {
    double directRad = k/180.*M_PI;

    for(int i=0; i<sizeof(fi0)/sizeof(double); ++i) {

      qDebug() << "--------------------------------";
      qDebug() << "direct=" << k << " {fi,lm}=" <<  fi0[i]/M_PI*180. << "," << lm0[i]/M_PI*180.;

      std::vector<SD::pointInArc> outputPoints;
      
      a.intersecare(fi0[i],lm0[i],directRad,&outputPoints);
      for(unsigned j=0; j<outputPoints.size(); ++j) {
        const SD::pointInArc &p = outputPoints[j];
        qDebug() << p.arcNum << ") " << p.geo.fi/M_PI*180 << p.geo.lm/M_PI*180;
      }
    }
  }
}

int main(int argc, char *argv[])
{
  test1();
  test2();
  return 0;
}
