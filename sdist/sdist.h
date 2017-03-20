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

#ifndef SDIST_H
#define SDIST_H
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <algorithm>

// Sphere Distance
namespace SD {

  class operations;

//класс точки декартовой, условные единицы
struct punctumXYZ {
  punctumXYZ(double x_, double y_, double z_)
    : x(x_), y(y_), z(z_) {}
  punctumXYZ()
      : x(0), y(0), z(0) {}
  double x,y,z;
  //
  double r2() {
    return x*x+y*y+z*z;
  }
};

//класс точки географической
struct punctumRad {
  punctumRad(double fi_, double lm_)
    : fi(fi_), lm(lm_){}
  punctumRad()
    : fi(0), lm(0){}
  //радианы
  double fi, lm;
};

//объединенная точка
struct punctum: punctumRad, punctumXYZ {
  punctum() {};
  punctum(const punctumRad &rad)
    : punctumRad(rad.fi, rad.lm) {}
  punctum(const punctumXYZ &xyz)
    : punctumXYZ(xyz.x, xyz.y, xyz.z) {}
  punctum(const punctumRad &rad, const punctumXYZ &xyz)
    : punctumRad(rad.fi, rad.lm) 
    , punctumXYZ(xyz.x, xyz.y, xyz.z) {}
};

//точка в контуре   
struct punctumOnAmbitus : punctum {
  //идентификатор контура
  int idAmbitus;
  //номер арки
  int idArcus;
  //расстояние в рад от базовой точки до точки контура
  double distRad;
  //
  punctumOnAmbitus &operator=(const punctumRad &p) {
    fi = p.fi;
    lm = p.lm;
    return *this;
  }
};

//класс плоскости проходящей через центр сферы
struct planum {
  //вектора нормали
  double a, b, c;

  //расчитать положение точки относительно плоскости
  double value(const punctumXYZ &inputPoint) const {
    return a*inputPoint.x+b*inputPoint.y+c*inputPoint.z;
  }
  //проверить лежит ли точка на плоскости
  bool contains(const punctumXYZ &inputPoint, const double &eps) {
    return a*inputPoint.x+b*inputPoint.y+c*inputPoint.z <= eps;
  }
};
  
//
class operations {
public:  
  //угол между направлением на север и дугой, построенной по двум точкам
  static double directionBy2Points(
    double fi1, double lm1, 
    double fi2, double lm2
  ) {
    
    if(fabs(lm1-lm2)<EPS) 
      return fi1<fi2 ? 0 : M_PI;
    double L = acos( cos(M_PI/2-fi1)*cos(M_PI/2-fi2) + 
                     sin(M_PI/2-fi1)*sin(M_PI/2-fi2)*cos(fabs(lm2-lm1))
      );

    double alfa = acos((cos(M_PI/2-fi2)-cos(M_PI/2-fi1)*cos(L))/(sin(M_PI/2-fi1)*sin(L)));
    return 
      lm1>lm2 ? 2*M_PI-alfa : alfa;
  }
  //угол между двумя плоскостями
  static double anglePlanumsRad(const planum &p1, const planum &p2) {
    double a = (p1.a*p2.a+p1.b*p2.b+p1.c*p2.c)/
      sqrt(p1.a*p1.a+p1.b*p1.b+p1.c*p1.c)/
        sqrt(p2.a*p2.a+p2.b*p2.b+p2.c*p2.c);
    return acos(a);
  }
  //перевод из радианов в градусы
  static double G(const double &rad) {
    return rad/M_PI*180;
  }
  //
  static punctumXYZ geo2xyz(const punctumRad &p) {
    punctumXYZ xyz;
    xyz.x = cos(p.fi)*cos(p.lm);
    xyz.y = cos(p.fi)*sin(p.lm);
    xyz.z = sin(p.fi);
    return xyz;
  }
  //
  static punctumRad xyz2geo(const punctumXYZ &p) {
    punctumRad g;
    g.fi = asin(p.z);
    double cos_fi = cos(g.fi);
    double cos_lm = p.x/cos_fi;
    double sin_lm = p.y/cos_fi;
    g.lm = atan2(sin_lm,cos_lm);
    return g;
  }
  //вычисление детерминанта 2
  //|a1 b1|
  //|a2 b2|
  //order a1,a2,b1,b2
  static double det2(double a1, double a2, 
    double b1, double b2) {
    return a1*b2-a2*b1;
  }
  //построение плоскости проходящей через центр сферы
  static planum planumCrossZero(
    const punctumXYZ &p1,
    const punctumXYZ &p2) {
    planum l;
    l.a =  det2(p1.y, p2.y, p1.z, p2.z);
    l.b = -det2(p1.x, p2.x, p1.z, p2.z);
    l.c =  det2(p1.x, p2.x, p1.y, p2.y);
    return l;
  }
  //пересечение 2 плоскостей проходящих через центр сферы 
  //и ограниченных сферой же
  //-1 - решения не существует
  //0 - плоскости совпадают
  static int intersecare(const planum &pm1, const planum &pm2, 
    punctumXYZ &xyz1, punctumXYZ &xyz2) {
    
    //может плоскости совпадают ?
    //проверяем тип коэффициентов
    int type1 = coefType(pm1);
    int type2 = coefType(pm2);

    if(type1==0 || type2==0)
      return -1;

    //|a1 b1|
    //|a2 b2|
    //order a1,a2,b1,b2
    double detAB = det2(pm1.a, pm2.a, pm1.b, pm2.b);

    //|a1 b1|
    //|a2 b2|
    //order a1,a2,b1,b2
    double detAC = det2(pm1.a, pm2.a, pm1.c, pm2.c);

    //|a1 b1|
    //|a2 b2|
    //order a1,a2,b1,b2
    double detBC = det2(pm1.b, pm2.b, pm1.c, pm2.c);

    //у второй плоскости все ненулевые коэффициенты
    if(type2==7) 
      return calc7(pm2,pm1,detAB,detAC,detBC,xyz1,xyz2);
    
    //у первой плоскости все ненулевые коэффициенты
    if(type1==7)
      return calc7(pm2,pm1,detAB,detAC,detBC,xyz1,xyz2);

    //плоскости совпадают
    if( (type1==1 && type2==1) ||
        (type1==2 && type2==2) ||
        (type1==3 && type2==3) ||
        ((type1==4 && type2==4) && fabs(detAB)<=EPS) ||
        ((type1==5 && type2==5) && fabs(detAC)<=EPS) ||
        ((type1==6 && type2==6) && fabs(detBC)<=EPS) 
    )
      return 0;

    if( (type1==1 && type2==3) ||
        (type1==1 && type2==4) ||
        (type1==3 && type2==4) ||
        (type2==1 && type1==3) ||
        (type2==1 && type1==4) ||
        (type2==3 && type1==4) ||
        ((type1==4 && type2==4) && fabs(detAB)>EPS)
    ) {
      xyz1.x = 0;   xyz2.x = 0; 
      xyz1.y = 0;   xyz2.y = 0;
      xyz1.z = 1.;  xyz2.z = -1.;
      return 2;
    }

    if( (type1==1 && type2==2) ||
        (type1==1 && type2==6) ||
        (type1==2 && type2==6) ||
        (type2==1 && type1==2) ||
        (type2==1 && type1==6) ||
        (type2==2 && type1==6) ||
        ((type1==6 && type2==6) && fabs(detBC)>EPS)
    ) {
      xyz1.x = 1.;  xyz2.x = -1.; 
      xyz1.y = 0;   xyz2.y = 0;
      xyz1.z = 0;   xyz2.z = 0;
      return 2;
    }

    if( (type1==2 && type2==3) ||
        (type1==2 && type2==5) ||
        (type1==3 && type2==5) ||
        (type2==2 && type1==3) ||
        (type2==2 && type1==5) ||
        (type2==3 && type1==5) ||
        ((type1==5 && type2==5) && fabs(detAC)>EPS)
    ) {
      xyz1.x = 0.;  xyz2.x = 0; 
      xyz1.y = 1;   xyz2.y = -1;
      xyz1.z = 0;   xyz2.z = 0;
      return 2;
    }

    if( (type1==1 && type2==5) ) {
      return calc15(pm1,pm2,xyz1,xyz2);
    }
    if( (type1==5 && type2==1) ) {
      return calc15(pm2,pm1,xyz1,xyz2);
    }

    if( (type1==2 && type2==4) ) {
      return calc24(pm1,pm2,xyz1,xyz2);
    }
    if( (type1==4 && type2==2) ) {
      return calc24(pm2,pm1,xyz1,xyz2);
    }

    if( (type1==3 && type2==6) ) {
      return calc36(pm1,pm2,xyz1,xyz2);
    }
    if( (type1==6 && type2==3) ) {
      return calc36(pm2,pm1,xyz1,xyz2);
    }

    if( (type1==4 && type2==5) ) {
      return calc45(pm1,pm2,xyz1,xyz2);
    }
    if( (type1==5 && type2==4) ) {
      return calc45(pm2,pm1,xyz1,xyz2);
    }

    if( (type1==4 && type2==6) ) {
      return calc46(pm1,pm2,xyz1,xyz2);
    }
    if( (type1==6 && type2==4) ) {
      return calc46(pm2,pm1,xyz1,xyz2);
    }

    if( (type1==5 && type2==6) ) {
      return calc56(pm1,pm2,xyz1,xyz2);
    }
    if( (type1==6 && type2==5) ) {
      return calc56(pm2,pm1,xyz1,xyz2);
    }
    return 0;
  }
  //
  static int calc7(const planum &pm1, const planum &pm2,
    double detAB, double detAC, double detBC,
    punctumXYZ &xyz1, punctumXYZ &xyz2) {

    if(fabs(detAB)<=EPS) {

      if(fabs(detAC)<=EPS) {
        //71
        return 0;
      }    
      //72         
      xyz1.z = 0;
      xyz1.y = +1./sqrt(1+(pm1.b/pm1.a)*(pm1.b/pm1.a));
      xyz1.x = -pm1.b/pm1.a * xyz1.y;

      xyz2.z = 0;
      xyz2.y = - xyz1.y;
      xyz2.x = - xyz1.x;
      return 2;
    } 
    
    //73
    double M = -detAC/detAB;
    double N = -(pm1.b/pm1.a)*M-(pm1.c/pm1.a);
    xyz1.z = +1 / sqrt(1+M*M+N*N);
    xyz1.x = xyz1.z * N;
    xyz1.y = xyz1.z * M;

    xyz2.z = -xyz1.z;
    xyz2.x = -xyz1.x;
    xyz2.y = -xyz1.y;
    return 2;
  }
  //
  static int calc15(const planum &pm1, const planum &pm2,
    punctumXYZ &xyz1, punctumXYZ &xyz2
  ) {    
    xyz1.z = +1./sqrt(1+(pm2.c/pm2.a)*(pm2.c/pm2.a));
    xyz1.x = -pm2.c/pm2.a * xyz1.z;
    xyz1.y = 0;

    xyz2.z = - xyz1.z;
    xyz2.x = - xyz1.x;
    xyz2.y = 0;
    return 2;
  }
  //
  static int calc24(const planum &pm1, const planum &pm2,
    punctumXYZ &xyz1, punctumXYZ &xyz2
  ) {    
    xyz1.y = +1./sqrt(1+(pm2.b/pm2.a)*(pm2.b/pm2.a));
    xyz1.x = -pm2.b/pm2.a * xyz1.y;
    xyz1.z = 0;

    xyz2.y = - xyz1.y;
    xyz2.x = - xyz1.x;
    xyz2.z = 0;
    return 2;
  }
  //
  static int calc36(const planum &pm1, const planum &pm2,
    punctumXYZ &xyz1, punctumXYZ &xyz2
  ) {    
    xyz1.z = +1./sqrt(1+(pm2.c/pm2.b)*(pm2.c/pm2.b));
    xyz1.y = -pm2.c/pm2.b * xyz1.z;
    xyz1.x = 0;

    xyz2.z = - xyz1.z;
    xyz2.y = - xyz1.y;
    xyz2.x = 0;
    return 2;
  }
  //
  static int calc45(const planum &pm1, const planum &pm2,
    punctumXYZ &xyz1, punctumXYZ &xyz2
  ) {    
    xyz1.x = +1./sqrt(1+(pm1.a/pm1.b)*(pm1.a/pm1.b)+(pm2.a/pm2.c)*(pm2.a/pm2.c));
    xyz1.y = -pm1.a/pm1.b * xyz1.x;
    xyz1.z = -pm2.a/pm2.c * xyz1.x;

    xyz2.x = - xyz1.x;
    xyz2.y = - xyz1.y;
    xyz2.z = - xyz1.z;
    return 2;
  }
  //
  static int calc46(const planum &pm1, const planum &pm2,
    punctumXYZ &xyz1, punctumXYZ &xyz2
  ) {    
    xyz1.y = +1./sqrt(1+(pm1.b/pm1.a)*(pm1.b/pm1.a)+(pm2.b/pm2.c)*(pm2.b/pm2.c));
    xyz1.x = -pm1.b/pm1.a * xyz1.y;
    xyz1.z = -pm2.b/pm2.c * xyz1.y;

    xyz2.y = - xyz1.y;
    xyz2.x = - xyz1.x;
    xyz2.z = - xyz1.z;
    return 2;
  }
  //
  static int calc56(const planum &pm1, const planum &pm2,
    punctumXYZ &xyz1, punctumXYZ &xyz2
  ) {    
    xyz1.z = +1./sqrt(1+(pm1.c/pm1.a)*(pm1.c/pm1.a)+(pm2.c/pm2.b)*(pm2.c/pm2.b));
    xyz1.x = -pm1.c/pm1.a * xyz1.z;
    xyz1.y = -pm2.c/pm2.b * xyz1.z;

    xyz2.z = - xyz1.z;
    xyz2.x = - xyz1.x;
    xyz2.y = - xyz1.y;
    return 2;
  }
  //-1 совпадают плоскости
  static int variant(const planum &pm1, const planum &pm2) {
    if( fabs(pm1.a-pm2.a)<EPS &&
        fabs(pm1.b-pm2.b)<EPS &&
        fabs(pm1.c-pm2.c)<EPS )
    return -1;

    //b1*y=0;
    if(fabs(pm1.a)<EPS && fabs(pm1.c)<EPS) {
    }

  }
  //тип коэффициентов плоскости
  static int coefType(const planum &pm) {
    //b*y=0;
    if(fabs(pm.a)<=EPS && fabs(pm.b)>EPS && fabs(pm.c)<=EPS) {
      return 1;
    }
    //c*z=0;
    if(fabs(pm.a)<=EPS && fabs(pm.b)<=EPS && fabs(pm.c)>EPS) {
      return 2;
    }
    //a*x=0;
    if(fabs(pm.a)> EPS && fabs(pm.b)<=EPS && fabs(pm.c)<=EPS) {
      return 3;
    }
    //a*x+b*y=0;
    if(fabs(pm.a)> EPS && fabs(pm.b)>EPS && fabs(pm.c)<=EPS) {
      return 4;
    }
    //a*x+c*z=0;
    if(fabs(pm.a)> EPS && fabs(pm.b)<=EPS && fabs(pm.c)>EPS) {
      return 5;
    }
    //b*y+c*z=0;
    if(fabs(pm.a)<= EPS && fabs(pm.b)>EPS && fabs(pm.c)>EPS) {
      return 6;
    }
    //0=0;
    if(fabs(pm.a)<=EPS && fabs(pm.b)<=EPS && fabs(pm.c)<=EPS) {
      return 0;
    }
    //все не нули
    return 7;
  }

  //расстояние в радианах между двумя точками
  static double distance(const punctumXYZ &xyz1, const punctumXYZ &xyz2) {
    double L = sqrt((xyz1.x-xyz2.x)*(xyz1.x-xyz2.x) +
                (xyz1.y-xyz2.y)*(xyz1.y-xyz2.y) +
                (xyz1.z-xyz2.z)*(xyz1.z-xyz2.z)
                );
    return 2*asin(L/2);
  };

  //нахождение точки по направлению на заданном расстоянии в рад
  static punctumRad toDirection(const punctumRad &g0,
    const double &directRad, const double &distanceRad);

  //проверить непрерывность пары точек
  static double continuity(const double &prev, const double &current) {
    double withShift = current + 2*M_PI;
    double delta1 = fabs(prev-current);
    double delta2 = fabs(prev-withShift);
    if(delta1<delta2)
      return current;
    else
      return withShift;
  }

  //строим перпендикулярные плоскости проходящие через m_pnt1, m_pnt2
  //разворачиваем их так чтобы при подстановке m_pnt1 в 2 и при подстановке m_pnt2 в 1 
  //выполнялось > 0
  static planum planumNormalCrossZero(const planum &plane, 
    const punctumXYZ &p1, const punctumXYZ &p2) {
    planum output;
    output.a =  det2(p1.y, plane.b, p1.z, plane.c);
    output.b = -det2(p1.x, plane.a, p1.z, plane.c);
    output.c =  det2(p1.x, plane.a, p1.y, plane.b);
    if(output.value(p2)<0) {      
      output.a *=-1; output.b *=-1; output.c *=-1;      
    }
    return output;
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
  static void operations::buildSector(
    const SD::punctumRad &g0,
    const double &Skm, 
    const double &startDirectRad,
    const double &stopDirectRad,
    double const &stepRad,
    const double &Re,
    std::vector<SD::punctumRad> &vSegment
  );

  //построение кругового сектора заданного размера
  //выполняется минимальная проверка на корректность входных данных
  //g0 - центр сектора
  //Skm - радиус сектора в км
  //startDirectRad - начальный направляющий угол сектора в рад (от направления на север по часовой стрелке)
  //stopDirectRad - конечный направляющий угол сектора в рад (от направления на север по часовой стрелке)
  //stepRad - шаг в рад, с которым строится ломаная
  //Re - радиус Земли в м
  //vSegment - замкнутая ломаная
  static void operations::buildSector2(
    const SD::punctumRad &g0,
    const double &Skm, 
    const double &startDirectRad,
    const double &stopDirectRad,
    double const &stepRad,
    const double &Re,
    std::vector<SD::punctumRad> &vSegment
  );

  static double EPS;
};

//класс арки
class arcus {
public:

  //
  arcus (const punctumXYZ &pnt1_, const punctumXYZ &pnt2_, const double &baseLm)
    : m_pnt1(pnt1_), m_pnt2(pnt2_), m_baseLm(baseLm) {
    m_planum = operations::planumCrossZero(m_pnt1,m_pnt2);

    //строим перпендикулярные плоскости проходящие через m_pnt1, m_pnt2
    //разворачиваем их так чтобы при подстановке m_pnt1 в 2 и при подстановке m_pnt2 в 1 
    //выполнялось > 0
    m_planum_p1 = operations::planumNormalCrossZero(m_planum, m_pnt1,m_pnt2);
    m_planum_p2 = operations::planumNormalCrossZero(m_planum, m_pnt2,m_pnt1);
  }

  //пересечение дуги с плоскостью проходящей через центр координат
  //на входе плоскость и точка от которой строилась плоскость
  bool intersecare(const planum &inputPlanum, 
    const punctumXYZ &inputPoint,
    punctumXYZ &outputXYZ) {
    
    punctumXYZ xyz1, xyz2;
    int ret = operations::intersecare(m_planum, inputPlanum, xyz1, xyz2);
    switch(ret) {
      case 0: {
        //плоскости совпадают - смотрим лежит ли исходная точка на дуге
        if(inArcus(inputPoint)) {
          outputXYZ = inputPoint;
          return true;
        }
        return false;
      }
      case -1:{
        return false;
      }
      case 1: {
        if(inArcus(xyz1)) {
          outputXYZ = xyz1;
          return true;
        } else
          return false;
      }//case 1
      case 2: {
        if(inArcus(xyz1)) {
          outputXYZ = xyz1;
          return true;
        } 
        if(inArcus(xyz2)) {
          outputXYZ = xyz2;
          return true;
        } 
        return false;
      }//case 2

    }//switch
    return false;
  }

  //
  const punctumXYZ *p1() const { return &m_pnt1; }
  const punctumXYZ *p2() const { return &m_pnt2; }
  const double baseLm() { return m_baseLm; }

  //определить лежит ли точка на дуге, ограниченной точками m_pnt1, m_pnt2
  bool inArcus(const punctumXYZ &xyz) const {
    return m_planum_p1.value(xyz) >= 0 && m_planum_p2.value(xyz) >= 0;
  }

private:
  //две точки в декартовых координатах
  const punctumXYZ m_pnt1, m_pnt2;
  planum m_planum; 
  //перепндикуляр через 1
  planum m_planum_p1; 
  //перепндикуляр через 2
  planum m_planum_p2; 
  //широта одного из концов арки для поддержки непрерывности
  double m_baseLm;
};

//класс контура
class ambitus {
public:
  /////////////////////////////////////////////////////////////////////
  ambitus(int id_) : idAmbitus(id_) {}
  /////////////////////////////////////////////////////////////////////
  //добавить точку в контур
  void append(const double &fi, const double &lm);
  //подготовить контур для расчетов
  //bClose - признак замыкания контура
  void prepare(const bool &bClose);
  /////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////
  //пересечение контура с большим кругом, проходящим через заданную точку в заданном направлении
  //направление расчитывается от направления на север по часовой стрелке
  //fiRad, lmRad - широта/долгота заданной точки в рад
  //fiRad=[-Пи/2,Пи/2], lmRad=[-Пи,Пи] - широта/долгота заданной точки в рад
  //directRad - заданное направление, радианы
  //outputPoints - массив точек пересечения
  //numberPunctumByDirection - номер точки в массиве outputPoints, самой первой от базовой
  //точки по заданному направлению
  //distanceRad - на сколько отступать в заданном направлении для построения большого круга
  //по умолчанию на 1 минуту
  void intersecare(
    const double &fiRad, const double &lmRad, 
    const double &directRad,    
    std::vector<punctumOnAmbitus> *outputPoints,
    int &numberPunctumByDirection,
    const double &distanceRad = M_PI/180./60.
  );
  /////////////////////////////////////////////////////////////////////
  //пересечение контура с дугой, заданной 2-мя точками
  //fiRad1, lmRad1 - широта/долгота 1-й точки в рад
  //fiRad2, lmRad2 - широта/долгота 1-й точки в рад
  //fiRad=[-Пи/2,Пи/2], lmRad=[-Пи,Пи] - широта/долгота заданной точки в рад
  //outputPoints - массив точек пересечения
  void intersecare(
    const double &fiRad1, const double &lmRad1, 
    const double &fiRad2, const double &lmRad2, 
    std::vector<punctumOnAmbitus> *outputPoints
  );
  /////////////////////////////////////////////////////////////////////
  //пересечение контура с плоскостью проходящей через центр координат
  //inputPlanum - плоскость пересечения
  //inputPoint - исходная точка, через которую проходит плоскость пересечения
  //outputPoints - массив точек пересечения
  void intersecare(
    const planum &inputPlanum, 
    const punctumXYZ &inputPoint,
    std::vector<punctumOnAmbitus> *outputPoints
  );
  /////////////////////////////////////////////////////////////////////
  //найти минимальное расстояние до контура - обход с шагом в градусы по кругу
  //медленная функция но надежная
  bool slowMinDistance(
    const punctumRad &inputPoint,
    punctumRad &outputPoint,
    double &distRad,
    const double &angleStepGrad = 1.
  );
  /////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////
  const std::vector<punctum> &puncta() const {
    return vPunctum;
  }

private:
  //идентификатор контура
  int idAmbitus;
  //массив точек
  std::vector<punctum> vPunctum;
  //арки из массива точек
  std::vector<arcus> vArcus;
};



};//

#endif // SDIST_H
