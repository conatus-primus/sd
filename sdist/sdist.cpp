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
//    Sources are located by https://github.com/conatus-primus/sd                       //
//                                                                                      //
//**************************************************************************************//

#include "sdist.h"


double SD::operations::EPS = 1e-10;

namespace gislab {
/****************************************************************************/
/*                          код предоставлен / начало                       */
/*      http://gis-lab.info/qa/sphere-geodesic-direct-problem.html          */  
/*                                                                          */
/****************************************************************************/
/* This file was automatically generated.  Do not edit! */
int SphereLinear(double pt1[],double pt2[],double dist13,double dist23,int clockwise,double pt3[]);
int SphereAngular(double pt1[],double pt2[],double azi13,double azi23,double pt3[]);
void SphereDirect(double pt1[],double azi,double dist,double pt2[]);
double CartToSpher(double x[],double y[]);
void Rotate(double x[],double a,int i);
void SpherToCart(double y[],double x[]);
void SphereInverse(double pt1[],double pt2[],double *azi,double *dist);
#define A_E 6371.0
#define Degrees(x) (x * 57.29577951308232)
#define Radians(x) (x / 57.29577951308232)

/****************************************************************************/
/*                          код предоставлен / окончание                    */
/*      http://gis-lab.info/qa/sphere-geodesic-direct-problem.html          */  
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                          код предоставлен / начало                       */
/*      http://gis-lab.info/qa/sphere-geodesic-direct-problem.html          */  
/*                                                                          */
/****************************************************************************/

#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <math.h>

void SphereInverse(double pt1[], double pt2[], double *azi, double *dist)
{
  double x[3], pt[2];

  SpherToCart(pt2, x);
  Rotate(x, pt1[1], 2);
  Rotate(x, M_PI_2 - pt1[0], 1);
  CartToSpher(x, pt);
  *azi = M_PI - pt[1];
  *dist = M_PI_2 - pt[0];

  return;
}

void SphereDirect(double pt1[], double azi, double dist, double pt2[])
{
  double pt[2], x[3];

  pt[0] = M_PI_2 - dist;
  pt[1] = M_PI - azi;
  SpherToCart(pt, x);
  Rotate(x, pt1[0] - M_PI_2, 1);
  Rotate(x, -pt1[1], 2);
  CartToSpher(x, pt2);

  return;
}

int SphereAngular(double pt1[], double pt2[], double azi13, double azi23,
		  double pt3[])
{
  double azi12, dist12, azi21, dist13;
  double cos_beta1, sin_beta1, cos_beta2, sin_beta2, cos_dist12, sin_dist12;

  SphereInverse(pt2, pt1, &azi21, &dist12);
  SphereInverse(pt1, pt2, &azi12, &dist12);
  cos_beta1 = cos(azi13 - azi12);
  sin_beta1 = sin(azi13 - azi12);
  cos_beta2 = cos(azi21 - azi23);
  sin_beta2 = sin(azi21 - azi23);
  cos_dist12 = cos(dist12);
  sin_dist12 = sin(dist12);

  if (sin_beta1 == 0. && sin_beta2 == 0.)     // Решение - любая точка
    return -1;				      //   на большом круге Q1-Q2.
  else if (sin_beta1 == 0.) {
    pt3[0] = pt2[0];			      // Решение - точка Q2.
    pt3[1] = pt2[1];
    return 0;
  } else if (sin_beta2 == 0.) {		      // Решение - точка Q1.
    pt3[0] = pt1[0];
    pt3[1] = pt1[1];
    return 0;
  } else if (sin_beta1 * sin_beta2 < 0.) {    // Лучи Q1-Q3 и Q2-Q3 направлены
    if (fabs(sin_beta1) >= fabs(sin_beta2)) { //   в разные полусферы.
      cos_beta2 = -cos_beta2;		      // Выберем ближайшее решение:
      sin_beta2 = -sin_beta2;		      //   развернём луч Q2-Q3 на 180°;
    } else {				      //     иначе
      cos_beta1 = -cos_beta1;		      //   развернём луч Q1-Q3 на 180°.
      sin_beta1 = -sin_beta1;
    }
  }
  dist13 = atan2(fabs(sin_beta2) * sin_dist12,
		 cos_beta2 * fabs(sin_beta1)
		 + fabs(sin_beta2) * cos_beta1 * cos_dist12);
  SphereDirect(pt1, azi13, dist13, pt3);

  return 0;
}

int SphereLinear(double pt1[], double pt2[], double dist13, double dist23,
		 int clockwise, double pt3[])
{
  double azi12, dist12, azi13;
  double cos_beta1;

  if (dist13 == 0.) {			      // Решение - точка Q1.
    pt3[0] = pt1[0];
    pt3[1] = pt1[1];
    return 0;
  } else if (dist23 == 0.) {		      // Решение - точка Q2.
    pt3[0] = pt2[0];
    pt3[1] = pt2[1];
    return 0;
  }

  SphereInverse(pt1, pt2, &azi12, &dist12);
  cos_beta1 = (cos(dist23) - cos(dist12) * cos(dist13))
    / (sin(dist12) * sin(dist13));
  if (fabs(cos_beta1) > 1.)		      // Решение не существует.
    return -1;
  azi13 = clockwise ? azi12 + acos(cos_beta1) : azi12 - acos(cos_beta1);
  SphereDirect(pt1, azi13, dist13, pt3);

  return 0;
}

void Rotate(double x[], double a, int i)
{
  double c, s, xj;
  int j, k;

  j = (i + 1) % 3;
  k = (i - 1) % 3;
  c = cos(a);
  s = sin(a);
  xj = x[j] * c + x[k] * s;
  x[k] = -x[j] * s + x[k] * c;
  x[j] = xj;

  return;
}

void SpherToCart(double y[], double x[])
{
  double p;

  p = cos(y[0]);
  x[2] = sin(y[0]);
  x[1] = p * sin(y[1]);
  x[0] = p * cos(y[1]);

  return;
}

double CartToSpher(double x[], double y[])
{
  double p;

  p = sqrt(x[0] * x[0] + x[1] * x[1]);
  y[1] = atan2(x[1], x[0]);
  y[0] = atan2(x[2], p);

  return sqrt(p * p + x[2] * x[2]);
}
/****************************************************************************/
/*                          код предоставлен / окончание                    */
/*      http://gis-lab.info/qa/sphere-geodesic-direct-problem.html          */  
/*                                                                          */
/****************************************************************************/
};

//нахождение точки по направлению на заданном расстоянии в рад
SD::punctumRad SD::operations::toDirection(const SD::punctumRad &g0,
  const double &directRad, const double &distanceRad) {
  //решаем прямую геодезическую задачу
  double pt1[] = {g0.fi, g0.lm};  
  double pt2[2];  
  gislab::SphereDirect(pt1,directRad,distanceRad,pt2);
  //проверяем на разрывность
  pt2[1] = operations::continuity(pt1[1],pt2[1]);
  return punctumRad(pt2[0],pt2[1]);
}


//класс контура
void SD::ambitus::append(const double &fi, const double &lm) {
  punctum point;
  point.fi = fi; 
  point.lm = lm;
  *(punctumXYZ*)(&point) = operations::geo2xyz(point);
  vPunctum.push_back(point);
}
//
//передается признак замыкания
void SD::ambitus::prepare(const bool &bClose) {
  const punctumRad &first = vPunctum.at(0);
  const punctumRad &last = vPunctum.at(vPunctum.size()-1);

  //замыкаем
  if(bClose && (fabs(first.fi-last.fi)>operations::EPS || fabs(first.lm-last.lm)>operations::EPS)) {
    vPunctum.push_back(vPunctum.at(0));
  }

  //формируем арки
  for(unsigned i=0; i<vPunctum.size()-1; i++) {
    arcus a(&vPunctum[i],&vPunctum[i+1]);
    vArcus.push_back(a);
  }
}
//пересечение контура с большим кругом проходящим через заданную точку в заданном направлении
//направление расчитывается от направления на север по часовой стрелке
//distanceRad - на сколько отступать в заданном направлении для построения большого круга
//по умолчанию на 1 минуту
void SD::ambitus::intersecare(
  const double &fiRad, const double &lmRad, 
  const double &directRad,    
  std::vector<punctumOnAmbitus> *outputPoints,
  const double &distanceRad
  ) {

  SD::punctumRad g0(fiRad,lmRad);
  SD::punctumRad g1 = SD::operations::toDirection(g0,directRad,distanceRad);

  //перейти в декартовы координаты
  SD::punctumXYZ xyz0 = SD::operations::geo2xyz(g0);
  SD::punctumXYZ xyz1 = SD::operations::geo2xyz(g1);

  //построить плоскость проходяющую через g0 и g1
  //построение плоскости проходящей через центр сферы
  SD::planum inputPlanum = SD::operations::planumCrossZero(xyz0, xyz1);
  punctumXYZ inputPoint;

  //пересечение контура с плоскостью проходящей через центр координат
  intersecare(inputPlanum, xyz0, outputPoints);
}
//пересечение контура с плоскостью проходящей через центр координат
void SD::ambitus::intersecare(const planum &inputPlanum, 
  const punctumXYZ &inputPoint,
  std::vector<punctumOnAmbitus> *outputPoints
  ) {
  bool bIntersect = false;
  for(unsigned  i=0; i<vArcus.size(); i++) {
    //пересечение дуги с плоскостью проходящей через центр координат
    punctumOnAmbitus output;
    if(!vArcus[i].intersecare(inputPlanum, inputPoint, output))
      continue;
    //
    output.distRad = operations::distance(inputPoint,output);
    output.idArcus = i;
    output.idAmbitus = idAmbitus;
    output = operations::xyz2geo(output);
    output.lm = operations::continuity(
      (vArcus[i].p1()->lm+vArcus[i].p2()->lm)/2., output.lm);
    outputPoints->push_back(output);
    bIntersect = true;
  }
}
//кол-во точек
int SD::ambitus::pointsCount() const {
  return vPunctum.size();
}
//кол-во арок
int SD::ambitus::arcusCount() const {
  return vArcus.size();
}
//получить одну точку по номеру
bool SD::ambitus::point(int number, punctum &p) const {
  if(number<vPunctum.size()) {
    p = vPunctum.at(number);
    return true;
  }
  return false;
}


//найти минимальное расстояние до контура 
bool SD::ambitus::slowMinDistance(
  const punctumRad &inputPoint,
  punctumRad &outputPoint,
  double &distRad,
  const double &angleStepGrad
) {
  std::vector<punctumOnAmbitus> outputPoints;
  for(double i=0; i<360; i+=angleStepGrad) {
    SD::ambitus::intersecare(
      inputPoint.fi, inputPoint.lm, 
      i/180.*M_PI,    
      &outputPoints
    );
  } 
  //расчитать расстояния и найти минимальное
  double distRadMin = DBL_MAX;  
  SD::punctumXYZ xyz0 = operations::geo2xyz(inputPoint);

  int number=0, numberMin=-1;
  for(const auto &point: outputPoints) {
    double distRad = operations::distance(point,xyz0);
    if(distRad<distRadMin) {
      distRadMin = distRad;
      numberMin = number;
    }
    number++;
  }
  if(numberMin!=-1) {
    outputPoint = outputPoints[numberMin];
    distRad = distRadMin;
    return true;
  }
  return false;
}

//сортировка точек по заданному направлению от базовой точки по возрастанию + 
//в противоположном направлении по возрастанию
void SD::operations::sort(const punctumXYZ &inputPointXYZ, 
  const std::vector<punctumOnAmbitus> &inputPuncta, 
  const double &directRad,    
  std::vector<punctumOnAmbitus> &outputPunctaProDirection, 
  std::vector<punctumOnAmbitus> &outputPunctaContraDirection
) {
  if(!inputPuncta.size())
    return;
  //находим самое большое расстояние
  int maxIndex = 0;
  double maxDistRad = inputPuncta.at(0).distRad;
  punctumOnAmbitus maxPunctum = inputPuncta.at(0);
  for(int i=1; i<inputPuncta.size(); i++) {
    if(maxDistRad<inputPuncta.at(i).distRad) {
      maxDistRad = inputPuncta.at(i).distRad;
      maxPunctum = inputPuncta.at(i);
      maxIndex = i;
    }
  }

  //смотрим точка лежит по заданному направлению или против
  punctumRad inputPointRad = operations::xyz2geo(inputPointXYZ);
  //нахождение точки по направлению на заданном расстоянии в рад
  SD::punctumRad punctumRadMaxPro = SD::operations::toDirection(
    inputPointRad, directRad, maxDistRad);

  SD::punctumRad punctumRadMaxContra = SD::operations::toDirection(
    inputPointRad, directRad+M_PI, maxDistRad);

  //смотрим совпадают координаты с исходными или сильно расходятся?
  SD::punctumXYZ punctumXYZMaxPro = operations::geo2xyz(punctumRadMaxPro);  
  double distPro = operations::distance(inputPuncta.at(maxIndex),punctumXYZMaxPro);

  SD::punctumXYZ punctumXYZMaxContra = operations::geo2xyz(punctumRadMaxContra);  
  double distContra = operations::distance(inputPuncta.at(maxIndex),punctumXYZMaxContra);

  std::vector<punctumOnAmbitus> *thisOutputPunctaProDirection = NULL;
  std::vector<punctumOnAmbitus> *thisOutputPunctaContraDirection = NULL;
  SD::punctumXYZ punctumXYZMax;

  if(distPro<distContra) {
    //максимально удаленная точка по направлению
    thisOutputPunctaProDirection    = &outputPunctaProDirection;
    thisOutputPunctaContraDirection = &outputPunctaContraDirection;
    punctumXYZMax = punctumXYZMaxPro;
  } else {
    //максимально удаленная точка против направления
    thisOutputPunctaProDirection    = &outputPunctaContraDirection;
    thisOutputPunctaContraDirection = &outputPunctaProDirection;
    punctumXYZMax = punctumXYZMaxContra;
  }

  //строим арку
  punctum p1(punctumXYZMax), p2(inputPointXYZ);
  arcus outputArc(&p1, &p2);

  //бежим по точкам и смотрим попадают они в арку (базовая точка-самая удаленная точка)
  for(const auto &p: inputPuncta) {
    if(outputArc.inArcus(p)) 
      thisOutputPunctaProDirection->push_back(p);
    else 
      thisOutputPunctaContraDirection->push_back(p);
  }
}

