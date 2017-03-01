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
  punctum(){}
  punctum(const punctumRad &rad)
    : punctumRad(rad.fi, rad.lm) {}
  punctum(const punctumXYZ &xyz)
    : punctumXYZ(xyz.x, xyz.y, xyz.z) {}
};

//точка на арке
struct punctumOnArcus : punctum {
  //номер арки
  int idArcus;
};

//точка в контуре
struct punctumOnAmbitus : punctumOnArcus {
  //идентификатор контура
  int idAmbitus;
  //расстояние в рад от базовой точки до точки контура
  double distRad;
  punctumOnAmbitus &operator=(const punctumRad &p) {
    fi = p.fi;
    lm = p.lm;
    return *this;
  }
};

//класс плоскости проходящей через центр сферы
struct planum {
  double a, b, c;
  //проверить лежит ли точка на плоскости
  bool contains(const punctumXYZ &inputPoint, const double &eps) {
    return a*inputPoint.x+b*inputPoint.y+c*inputPoint.z <= eps;
  }
};

//
class operations {
public:
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

  //проверка на попадание точки внтурь или на концы отрезка
  static bool pointInSegmentXYZ(const punctumXYZ &xyz1, const punctumXYZ &xyz2, const punctumXYZ &inputXYZ) {
    
    double res;

    res = (inputXYZ.x-xyz1.x)*(inputXYZ.x-xyz2.x);    
    if(fabs(res)>EPS) {
      if(res>0)
        return false;
    }
    
    res = (inputXYZ.y-xyz1.y)*(inputXYZ.y-xyz2.y);
    if(fabs(res)>EPS) {
      if(res>0)
        return false;
    }

    res = (inputXYZ.z-xyz1.z)*(inputXYZ.z-xyz2.z);
    if(fabs(res)>EPS) {
      if(res>0)
        return false;
    }

    return true;
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
  //сортировка точек по заданному направлению от базовой точки по возрастанию + 
  //в противоположном направлении по возрастанию
  static void sort(const punctumXYZ &inputPoint, 
    const std::vector<punctumOnAmbitus> &inputPuncta, 
    const double &directRad,    
    std::vector<punctumOnAmbitus> &outputPunctaProDirection, 
    std::vector<punctumOnAmbitus> &outputPunctaContraDirection
  );
  //сортировка точек по заданному направлению от базовой точки по возрастанию + 
  //в противоположном направлении по возрастанию
  static void sort(
    const double &fiRad, 
    const double &lmRad,     
    const std::vector<punctumOnAmbitus> &inputPuncta, 
    const double &directRad,    
    std::vector<punctumOnAmbitus> &outputPunctaProDirection, 
    std::vector<punctumOnAmbitus> &outputPunctaContraDirection
  ) {
    punctumXYZ inputPoint = operations::geo2xyz(punctumRad(fiRad,lmRad));
    sort(inputPoint, 
        inputPuncta, 
        directRad,    
        outputPunctaProDirection, 
        outputPunctaContraDirection
    );
  }

  static double EPS;
};

//класс арки
class arcus {
public:

  //
  arcus (const punctum *pnt1_, const punctum *pnt2_)
    : m_pnt1(pnt1_), m_pnt2(pnt2_) {
    m_planum = operations::planumCrossZero(*m_pnt1,*m_pnt2);
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
  const punctum *p1() const { return m_pnt1; }
  const punctum *p2() const { return m_pnt2; }

  //определить лежит ли точка на дуге
  bool inArcus(const punctumXYZ &xyz) {
    //смотрим на проекционные координаты и считаем простое вхождение точки в отрезок
    //считаем что точка лежит на плоскости и на большой дуге!
    return operations::pointInSegmentXYZ(*m_pnt1, *m_pnt2, xyz);
  }

private:
  const punctum *m_pnt1, *m_pnt2;
  planum m_planum; 
};

//класс контура
class ambitus {
public:
  /////////////////////////////////////////////////////////////////////
  ambitus(int id_) : idAmbitus(id_) {}
  //
  void append(const double &fi, const double &lm);
  //
  //передается признак замыкания
  void prepare(const bool &bClose);
  /////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////
  //пересечение контура с большим кругом проходящим через заданную точку в заданном направлении
  //направление расчитывается от направления на север по часовой стрелке
  //distanceRad - на сколько отступать в заданном направлении для построения большого круга
  //по умолчанию на 1 минуту
  void intersecare(
    const double &fiRad, const double &lmRad, 
    const double &directRad,    
    std::vector<punctumOnAmbitus> *outputPoints,
    const double &distanceRad = M_PI/180./60.
   );
  //пересечение контура с плоскостью проходящей через центр координат
  void intersecare(const planum &inputPlanum, 
    const punctumXYZ &inputPoint,
    std::vector<punctumOnAmbitus> *outputPoints
   );
  /////////////////////////////////////////////////////////////////////

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
  //кол-во точек в контуре
  int pointsCount() const;
  //кол-во арок
  int arcusCount() const;
  //получить одну точку по номеру
  bool point(int number, punctum &p) const;
  /////////////////////////////////////////////////////////////////////

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
