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
#include <iostream>

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
    arcus a(vPunctum[i],vPunctum[i+1],vPunctum[i].lm);
    vArcus.push_back(a);
  }
}

// Return whether first element is greater than the second
bool less ( const SD::punctumOnAmbitus &elem1, const SD::punctumOnAmbitus &elem2 )
{
  return elem1.distRad < elem2.distRad;
}
/////////////////////////////////////////////////////////////////////
//пересечение контура с дугой, заданной 2-мя точками
//fiRad1, lmRad1 - широта/долгота 1-й точки в рад
//fiRad2, lmRad2 - широта/долгота 1-й точки в рад
//fiRad=[-Пи/2,Пи/2], lmRad=[-Пи,Пи] - широта/долгота заданной точки в рад
//outputPoints - массив точек пересечения
void SD::ambitus::intersecare(
  const double &fiRad1, const double &lmRad1, 
  const double &fiRad2, const double &lmRad2, 
  std::vector<punctumOnAmbitus> *outputPoints
  ) {
  if(!outputPoints)
    return;
  //строим плоскость через 2 точки
  punctumRad g1(fiRad1,lmRad1), g2(fiRad2,lmRad2);
  punctumXYZ pnt1 = operations::geo2xyz(g1), pnt2 = operations::geo2xyz(g2);

  //построение плоскости проходящей через центр сферы
  planum inputPlanum = operations::planumCrossZero(pnt1, pnt2);

  //находим точки пересечения
  intersecare(inputPlanum, pnt1, outputPoints);

  //находим точку ближайшую к начальной по дуге
  arcus a(pnt1,pnt2,g1.lm);
  std::vector<punctumOnAmbitus>::iterator it = outputPoints->begin();
  for(; it != outputPoints->end(); ) {
    if(a.inArcus(*it)) {
      it++;
      continue;
    }
    it = outputPoints->erase(it);
  }
  
  if(!outputPoints->size())
    return;

  //еще раз пробежимся по массиву и найдем ближайшую точку к 1-й
  std::sort(outputPoints->begin(),outputPoints->end(),less);
}

//пересечение контура с большим кругом проходящим через заданную точку в заданном направлении
//направление расчитывается от направления на север по часовой стрелке
//distanceRad - на сколько отступать в заданном направлении для построения большого круга
//по умолчанию на 1 минуту
void SD::ambitus::intersecare(
  const double &fiRad, const double &lmRad, 
  const double &directRad,    
  std::vector<punctumOnAmbitus> *outputPoints,
  int &numberPunctumByDirection,
  const double &distanceRad
  ) {
  if(!outputPoints)
    return;

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
  if(!outputPoints->size())
    return;

  //бежим в заданном направлении строим небольшие дуги и ищем вхождение в дугу первой
  //же точки, чтобы отследить самое первое пересечение
  double lmPrev = lmRad;
  punctum pm0(g0,xyz0);
  numberPunctumByDirection = -1;

  for(double i=0.5; i<180; i+=0.5) {
    double distRad = DBL_MAX;
    SD::punctumRad current = SD::operations::toDirection(
      g0, directRad,i/180.*M_PI);
    current.lm = SD::operations::continuity(lmPrev,current.lm);
    punctum pm1(current,operations::geo2xyz(current));
    arcus a(pm0,pm1,pm0.lm);
    for(unsigned k=0; k<outputPoints->size(); ++k) {
      const punctumOnAmbitus &pnctm = outputPoints->at(k);
      if(!a.inArcus(pnctm)) 
        continue;
      if(distRad<pnctm.distRad)
        continue;
      distRad = pnctm.distRad;
      numberPunctumByDirection = k;
    }
    if(numberPunctumByDirection!=-1) {
      break;
    }
    lmPrev = current.lm;
    pm0 = pm1;
  }
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
    output.lm = operations::continuity(vArcus[i].baseLm(), output.lm);
    outputPoints->push_back(output);
    bIntersect = true;
  }
}
//найти минимальное расстояние до контура 
bool SD::ambitus::slowMinDistance(
  const punctumRad &inputPoint,
  punctumRad &outputPoint,
  double &distRad,
  const double &angleStepGrad
) {
  std::vector<punctumOnAmbitus> outputPoints;
  int numberPunctumByDirection;
  for(double i=0; i<360; i+=angleStepGrad) {
    SD::ambitus::intersecare(
      inputPoint.fi, inputPoint.lm, 
      i/180.*M_PI,    
      &outputPoints,
      numberPunctumByDirection
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

//алгоритм построения выпуклой оболочки Джарвиса
namespace jarvise {

struct vv {
  double x, y;
  int index;
  vv(int index_, double x_, double y_) 
    : index(index_), x(x_), y(y_) {}
};

//векторное произведение
double vect(const vv &a1, const vv &a2, const vv &b1, const vv &b2) {
  return (a2.x - a1.x)*(b2.y-b1.y)-(b2.x-b1.x)*(a2.y-a1.y);
}
//квадратное расстояние
double dist2(const vv &a1, const vv &a2) {
  return sqrt(a2.x-a1.x)+sqrt(a2.y-a1.y);
}
//построение минимальной выпуклой оболочки
void solve(std::vector<vv> &a, int &k, std::vector<vv> &b) {
  int i, j, m;

  for( i= 0; i < a.size(); i++) {
    b.push_back(a[i]);
  }

  //ищем правую нижнюю точку
  m=0;
  for( i= 1; i < a.size(); i++) {
    if( a[i].y < a[m].y )
      m = i;
    else {
      if( (a[i].y == a[m].y) && (a[i].x > a[m].x) )
        m=i;
    }
  }

  //запишем ее в массив b и переставим на первое место в массиве a
  b[0] = a[m];   
  a[m] = a[0];   
  a[0] = b[0];
  k = 0;
  int min = 1;

  do  {
    //ищем очередную вершину оболочки
    for( j=1; j<a.size(); j++) {
      if ((vect(b[k],a[min],b[k],a[j])<0) ||
         (  (vect(b[k],a[min],b[k],a[j])<=SD::operations::EPS) &&
            (dist2(b[k],a[min])< dist2(b[k],a[j]))
         ))
        min=j;
    }
    k++;
    if(k==a.size())
      break;
    //записана очередная вершина
    b[k]=a[min];
    min=0;
  }
  while(!((b[k].x == b[0].x) && (b[k].y == b[0].y))); //пока ломаная не замкнется
}

};

//построение кругового сектора заданного размера
//выполняется минимальная проверка на корректность входных данных
//g0 - центр сектора
//Skm - радиус сектора в км
//startDirectRad - начальный направляющий угол сектора в рад (от направления на север по часовой стрелке)
//stopDirectRad - конечный направляющий угол сектора в рад (от направления на север по часовой стрелке)
//stepRad - шаг в рад, с которым строится ломаная
//Re - радиус Земли в м
//vSegment - замкнутая ломаная
void SD::operations::buildSector2(
  const SD::punctumRad &g0,
  const double &Skm, 
  const double &startDirectRad,
  const double &stopDirectRad,
  double const &stepRad,
  const double &Re,
  std::vector<SD::punctumRad> &vSegment
) {

  vSegment.clear();

  //минимальная проверка на корректность
  if((Skm==0) || (stepRad<operations::EPS) ||
    (startDirectRad==stopDirectRad)
  )
    return;

  //добавили центральную точку
  vSegment.push_back(g0);

  //строим длину сегмента
  double dbSegment = Skm*1000./Re;
    
  //добавили левый радиус
  for(double d=0; d<dbSegment; d+=stepRad) {
    SD::punctumRad currentPunctum = SD::operations::toDirection(
      g0, startDirectRad,d);
    vSegment.push_back(currentPunctum);
  }
  
  //точка начала дуги
  SD::punctumRad lPunctum = SD::operations::toDirection(
    g0, startDirectRad, dbSegment);

  //добавили последнюю точку левого радиуса и первую точку дуги
  vSegment.push_back(lPunctum);

  //точка конца дуги
  SD::punctumRad rPunctum = SD::operations::toDirection(
    g0, stopDirectRad, dbSegment);

  //расчитываем угол
  double angleSegment = SD::operations::directionBy2Points(
    lPunctum.fi, lPunctum.lm,
    rPunctum.fi, rPunctum.lm
  );
  
  //расчитываем длину арки
  double segmentArcLen = SD::operations::distance(
    SD::operations::geo2xyz(lPunctum),
    SD::operations::geo2xyz(rPunctum)
  );
  
  //дуга арки
  int N = 20;
  double stepRad1 = stepRad;
  if(segmentArcLen/stepRad<N)
    stepRad1 = segmentArcLen/N;
  for(double d=0; d<segmentArcLen; d+=stepRad1) {
    SD::punctumRad currentPunctum = SD::operations::toDirection(
      lPunctum, angleSegment,d);
    vSegment.push_back(currentPunctum);
  }
  
  //добавили последнюю точку дуги
  vSegment.push_back(rPunctum);

  //считаем правый радиус
  for(double d=dbSegment; d>=0; d-=stepRad) {
    SD::punctumRad currentPunctum = SD::operations::toDirection(
      g0, stopDirectRad,d);
    vSegment.push_back(currentPunctum);
  }

  //замыкаем центральной точкой
  vSegment.push_back(g0);
}

//построение кругового сектора заданного размера
//выполняется минимальная проверка на корректность входных данных
//g0 - центр сектора
//Skm - радиус сектора в км
//startDirectRad - начальный направляющий угол сектора в рад (от направления на север по часовой стрелке)
//stopDirectRad - конечный направляющий угол сектора в рад (от направления на север по часовой стрелке)
//stepRad - шаг в рад, с которым строится ломаная
//Re - радиус Земли в м
//vSegment - замкнутая ломаная
void SD::operations::buildSector(
  const SD::punctumRad &g0,
  const double &Skm, 
  const double &startDirectRad,
  const double &stopDirectRad,
  double const &stepRad,
  const double &Re,
  std::vector<SD::punctumRad> &vSegment
) {

  vSegment.clear();

  //минимальная проверка на корректность
  if((Skm==0) || (stepRad<operations::EPS) ||
    (startDirectRad==stopDirectRad)
  )
    return;

  int N = 20;

  //строим длину сегмента
  double dbSegment = Skm*1000./Re;

  //точка начала дуги
  SD::punctumRad lPunctum = SD::operations::toDirection(
    g0, startDirectRad, dbSegment);

  //точка конца дуги
  SD::punctumRad rPunctum = SD::operations::toDirection(
    g0, stopDirectRad, dbSegment);

  //точка середины дуги
  SD::punctumRad mPunctum = SD::operations::toDirection(
    g0, (stopDirectRad+startDirectRad)/2, dbSegment);

  //добавили центральную точку
  vSegment.push_back(g0);
    
  //добавили левый радиус
  for(double d=0; d<dbSegment; d+=stepRad) {
    SD::punctumRad currentPunctum = SD::operations::toDirection(
      g0, startDirectRad,d);
    vSegment.push_back(currentPunctum);
  }
  
  //добавили последнюю точку левого радиуса и первую точку дуги
  vSegment.push_back(lPunctum);

  //расчитываем угол
  double angleSegment1 = SD::operations::directionBy2Points(
    lPunctum.fi, lPunctum.lm,
    mPunctum.fi, mPunctum.lm
  );
  //расчитываем длину арки
  double segmentArcLen1 = SD::operations::distance(
    SD::operations::geo2xyz(lPunctum),
    SD::operations::geo2xyz(mPunctum)
  );

  //дуга арки
  double stepRad1 = stepRad;
  if(segmentArcLen1/stepRad<N)
    stepRad1 = segmentArcLen1/N;
  for(double d=0; d<segmentArcLen1; d+=stepRad1) {
    SD::punctumRad currentPunctum = SD::operations::toDirection(
      lPunctum, angleSegment1,d);
    vSegment.push_back(currentPunctum);
  }
  //добавили последнюю точку дуги
  vSegment.push_back(mPunctum);

  //расчитываем угол
  double angleSegment2 = SD::operations::directionBy2Points(
    mPunctum.fi, mPunctum.lm,
    rPunctum.fi, rPunctum.lm
  );
  //расчитываем длину арки
  double segmentArcLen2 = SD::operations::distance(
    SD::operations::geo2xyz(mPunctum),
    SD::operations::geo2xyz(rPunctum)
  );
    
  //дуга арки
  double stepRad2 = stepRad;
  if(segmentArcLen2/stepRad<N)
    stepRad2 = segmentArcLen2/N;
  for(double d=0; d<segmentArcLen2; d+=stepRad1) {
    SD::punctumRad currentPunctum = SD::operations::toDirection(
      mPunctum, angleSegment2,d);
    vSegment.push_back(currentPunctum);
  }  
  //добавили последнюю точку дуги
  vSegment.push_back(rPunctum);

  //считаем правый радиус
  for(double d=dbSegment; d>=0; d-=stepRad) {
    SD::punctumRad currentPunctum = SD::operations::toDirection(
      g0, stopDirectRad,d);
    vSegment.push_back(currentPunctum);
  }

  //замыкаем центральной точкой
  vSegment.push_back(g0);
}

