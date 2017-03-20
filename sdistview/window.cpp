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

#include <QtWidgets>
#include <QMenuBar>
#include <QMenu>
#include "C:/OSGeo4W/include/proj_api.h"
#include "renderarea.h"
#include "window.h"
#include "../sdist/sdist.h"


//радиус земли в м
#define Re (6371*1000)

//! [0]
const int IdRole = Qt::UserRole;
//! [0]

projPJ projSource, projTarget;

bool fl2pl(double fi, double lm,
  double &xM, double &yM) {
  xM = Re*lm; 
  yM = Re*sin(fi);
  
  //double lon=37.617778/180*M_PI;
  //double lat=55.751667/180*M_PI;
  yM = lm;
  xM = fi;
  double lon = lm;
  double lat = fi;
  int ret = pj_transform( projSource, projTarget, 1, 1,
                    &lon, &lat, NULL );
  if(ret<0) {
    char *err = pj_strerrno(ret);
  }
  xM = lon; yM = lat;
  return true;
}
bool pl2fl(double xM, double yM,
  double &fi, double &lm) {

  double lon = xM; 
  double lat = yM;

  int ret = pj_transform( projTarget, projSource, 1, 1,
                    &lon, &lat, NULL );
  if(ret<0) {
    char *err = pj_strerrno(ret);
  }
  lm = lon;
  fi = lat;

  //lm = xM/Re; 
  //fi = asin(yM/Re);
  //lm = xM; 
  //fi = yM;
  return true;
}

void f() {

  for(int a=1; a<=9; a++)
    for(int b=0; b<=9; b++)
      for(int c=0; c<=9; c+=2)
        for(int d=1; d<=9; d++)
          for(int e=0; e<=9; e++)
            for(int f=0; f<=9; f++)
              for(int g=0; g<=9; g++) {
                uint64_t x = a*100000 + b*10000+2*1000+3*100+40+c;
                uint64_t y = d*10000+e*100+f*10+g;
                if(y*72==x) {
                  //qDebug() << x << y;
                  break;
                }
              }

}

void f1() {
for(int a=1; a<=9; a++)
    for(int b=1; b<=9; b++)
      for(int c=1; c<=9; c++)
        for(int d=1; d<=9; d++)
          for(int e=1; e<=9; e++)
            for(int f=1; f<=9; f++)
              for(int g=1; g<=9; g++) 
                for(int x=1; x<=9; x++) 
                  for(int y=1; y<=9; y++) 
              {
                int m[9];
                memset(m,0,sizeof(m));
                if(a==b || a==c || a==d || a==e || a==f || a==g || a==x || a==y ||
                           b==c || b==d || b==e || b==f || b==g || b==x || b==y ||
                                   c==d || c==e || c==f || c==g || c==x || c==y ||
                                           d==e || d==f || d==g || d==x || d==y ||
                                                   e==f || e==g || e==x || e==y ||
                                                           f==g || f==x || f==y ||
                                                                   g==x || g==y ||
                                                                           x==y )
                  continue;

                int sum1 = a+b+c;
                int sum2 = d+e+f;
                int sum3 = g+x+y;
                int sum4 = a+d+g;
                int sum5 = b+e+x;
                int sum6 = c+f+y;
                int sum7 = g+e+c;
                int sum8 = a+e+y;
                if(
                  sum1==sum2 &&
                  sum1==sum3 &&
                  sum1==sum4 &&
                  sum1==sum5 &&
                  sum1==sum6 &&
                  sum1==sum7 &&
                  sum1==sum8 ) {

                    qDebug() << "-------------------------------";
                    qDebug() << a << b << c;
                    qDebug() << d << e << f;
                    qDebug() << g << x << y;
                }
              }
}

//
void Window::loadDialogFile() {

 m_stAmbitusFile = QFileDialog::getOpenFileName(this,
     tr("Open text"), "", tr("text Files (*.txt)"));

 
 loadFile();
 resizeData(); 
 update();
}

//
void Window::loadFile() {

  qDeleteAll(m_lAmbitus.begin(), m_lAmbitus.end());
  m_lAmbitus.clear();

  QFile file(m_stAmbitusFile);
  if(!file.open(QIODevice::ReadOnly|QIODevice::Text))
    return;

  QTextStream in(&file);
  ambitusPaint *pA  = NULL;
  int index = 0;
  while (!in.atEnd()) {

    QString line = in.readLine().trimmed();
    if(line.isEmpty())
      continue;

    QStringList l = line.split(" ");
    if(l.size()!=2) {
      if(pA) {
        m_lAmbitus.push_back(pA);
        pA = NULL;
      }
      continue;
    }

    if(l.size()==2) {
      //начало контура
      if(!pA) 
        pA = new ambitusPaint(index++);
      double fiCurrent = l.first().toDouble();
      double lmCurrent = l.last().toDouble();
        pA->a.append(fiCurrent, lmCurrent);
      if(pA->a.puncta().size()==1) {
        pA->fiMax = pA->fiMin = fiCurrent;
        pA->lmMax = pA->lmMin = lmCurrent;
      } else {
        pA->fiMax = qMax(pA->fiMax,fiCurrent);
        pA->fiMin = qMin(pA->fiMin,fiCurrent);
        pA->lmMax = qMax(pA->lmMax,lmCurrent);
        pA->lmMin = qMin(pA->lmMin,lmCurrent);
      }
    }

  }

  if(pA) {
    m_lAmbitus.push_back(pA);
  }
  if(!m_lAmbitus.size())
    return;

  m_fiMin = m_lAmbitus.at(0)->fiMin;
  m_fiMax = m_lAmbitus.at(0)->fiMax;
  m_lmMin = m_lAmbitus.at(0)->lmMin;
  m_lmMax = m_lAmbitus.at(0)->lmMax;

  m_lAmbitus[0]->prepare();
  for(int i=1; i<m_lAmbitus.size(); ++i) {
    m_fiMin = qMin(m_lAmbitus.at(i)->fiMin,m_fiMin);
    m_fiMax = qMax(m_lAmbitus.at(i)->fiMax,m_fiMax);
    m_lmMin = qMin(m_lAmbitus.at(i)->lmMin,m_lmMin);
    m_lmMax = qMax(m_lAmbitus.at(i)->lmMax,m_lmMax);
    m_lAmbitus[i]->prepare();
  }

  double delta = m_fiMax-m_fiMin;
  m_fiMin -= delta*0.3;
  m_fiMax += delta*0.3;
  if(m_fiMax>M_PI/2)
    m_fiMax = M_PI/2-0.01;
  
  delta = m_lmMax-m_lmMin;
  m_lmMin -= delta*0.1;
  m_lmMax += delta*0.1;

  /*
  m_fiMax = 89./180*M_PI;
  m_fiMin = -m_fiMax;
  m_lmMax = 2*M_PI;
  m_lmMin = -m_lmMax;
  */
  qDebug() << m_lAmbitus.size();
  //m_lmMax = M_PI-0.01;
}

//! [1]
Window::Window()
{ 
  projSource = pj_init_plus("+proj=longlat +ellps=WGS84"); 
  projTarget = pj_init_plus("+proj=longlat +ellps=WGS84"); 
  //projTarget = pj_init_plus(" +proj=longlat"); 
  //projTarget = pj_init_plus("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
  //projTarget = pj_init_plus("+proj=lcc +lat_1=35.25 +lat_2=36 +lat_0=34 +units=m"); 
  //projTarget = pj_init_plus("+proj=gnom +lat_0=90 +lon_0=0 +x_0=6300000 +y_0=6300000 +ellps=WGS84 +datum=WGS84 +units=m +no_defs");
  //projTarget = pj_init_plus("+proj=eqdc +lat_1=55 +lat_2=60");
  //projTarget = pj_init_plus("+proj=fouc"); 

  m_strobGrad = 10;
  m_Lkm=1700;
  m_Skm=800;
  m_fiCurrent = 55./180.*M_PI;
  m_lmCurrent = 37./180.*M_PI;
  m_directRad = 5.;
  m_numberPunctumByDirection = -1;

  setMouseTracking(true);
  setFocusPolicy(Qt::StrongFocus);
  m_bDeleteKeyPressed = false;

  update();
 
  QMenuBar* mnuBar = menuBar();
  QMenu*   pmnu   = new QMenu("&Menu");

  pmnu->addAction("Load", this, SLOT(loadDialogFile()));

  pmnu->addSeparator();

  QAction* pCheckableAction = pmnu->addAction("Contours with lines");
  pCheckableAction->setCheckable(true);
  pCheckableAction->setChecked(true);

  QMenu* pmnuSubMenu = new QMenu("&SubMenu", pmnu);
  pmnu->addMenu(pmnuSubMenu);
  pmnuSubMenu->addAction("&Test");

  QAction* pDisabledAction = pmnu->addAction("&DisabledItem");
  pDisabledAction->setEnabled(false);

  pmnu->addSeparator();

  pmnu->addAction("&Exit", qApp, SLOT(quit()));

  mnuBar->addMenu(pmnu);
  mnuBar->show();
}
//! [10]

Window::~Window() {
}

void Window::recalc() {
  m_vOutputPoints.clear();
  m_numberPunctumByDirection = -1;
  for(int i=0; i<m_lAmbitus.size(); ++i) {
    m_lAmbitus[i]->a.intersecare(
      m_fiCurrent, m_lmCurrent,     
      m_directRad/180.*M_PI,
      &m_vOutputPoints,
      m_numberPunctumByDirection
    );
  }
}

void Window::mouseDoubleClickEvent(QMouseEvent * event) {    
  xy2fl(event->x(), event->y(), m_fiCurrent, m_lmCurrent); 
  recalc();
  /*
  QTime tm;
  tm.start();
  m_vMinDistatncePoint.clear();
  double distRad;
  QMap<double,SD::punctumRad> dists;

  for(int i=0; i<m_lAmbitus.size(); ++i) {
    SD::punctumRad outPoint;
    if(!m_lAmbitus[i]->a.slowMinDistance(
      SD::punctumRad(m_fiCurrent, m_lmCurrent),
      outPoint,
      distRad
    ))
      continue;
    dists[distRad] = outPoint;
  }
  if(dists.size()) {
    m_vMinDistatncePoint << dists.first();
  }
  qDebug() << tm.elapsed();
*/  
  repaint();
}

void Window::mousePressEvent(QMouseEvent * event) {
}

void Window::keyReleaseEvent(QKeyEvent * event) {
  if(event->key()==Qt::Key_Delete || event->key()==Qt::Key_D)
    m_bDeleteKeyPressed = false;
  QWidget::keyReleaseEvent(event);
}

void Window::keyPressEvent(QKeyEvent * event) {

  if(event->key()==Qt::Key_X) {
    m_strobGrad += 1.;    
    if(m_strobGrad > 45)
      m_strobGrad = 1;
  }
  if(event->key()==Qt::Key_Z) {
    m_strobGrad -= 1.;    
    if(m_strobGrad <1)
      m_strobGrad = 45;
  }

  if(event->key()==Qt::Key_Up) {
    m_Lkm += 50.;    
    if(m_Lkm > 5000)
      m_Lkm = 50;
  }
  if(event->key()==Qt::Key_Down) {
    m_Lkm -= 50.;    
    if(m_Lkm < 0)
      m_Lkm = 5000;
  }

  if(event->key()==Qt::Key_Left) 
    m_directRad -= 1.;    
  if(event->key()==Qt::Key_Right)
    m_directRad += 1.;
  if(m_directRad>=360)
    m_directRad-=360;
  if(m_directRad<0)
    m_directRad+=360;

  if( event->key()==Qt::Key_Left || event->key()==Qt::Key_Right ||
      event->key()==Qt::Key_Up || event->key()==Qt::Key_Down ||
      event->key()==Qt::Key_X || event->key()==Qt::Key_Z 
  ) { 
    m_vOutputPoints.clear();
    m_numberPunctumByDirection = -1;
    for(int i=0; i<m_lAmbitus.size(); ++i) {
      m_lAmbitus[i]->a.intersecare(
        m_fiCurrent, m_lmCurrent,     
        m_directRad/180.*M_PI,
        &m_vOutputPoints,
        m_numberPunctumByDirection
      );
    }
    repaint();
  }

  if(event->key()==Qt::Key_Delete || event->key()==Qt::Key_D)
    m_bDeleteKeyPressed = true;
  QWidget::keyReleaseEvent(event);
}

void Window::mouseMoveEvent(QMouseEvent * event) {
  double fiCurrent, lmCurrent;
  xy2fl(event->x(), event->y(), fiCurrent, lmCurrent); 
  statusBar()->showMessage(tr("fi=%1 lm=%2").arg(fiCurrent/M_PI*180).arg(lmCurrent/M_PI*180));
}

QPolygon m_p;
void Window::paintEvent(QPaintEvent * event) {

  if(!m_lAmbitus.size())
    return;
  QPainter painter(this);

  //////////////////////////////////////////////////////////////////
  //нарисовали контуры
  QPen penContour(Qt::blue);
  penContour.setWidth(2);
  painter.setPen(penContour);
  //painter.setPen(QPen(QColor(Qt::blue)));
  int index = 0;
  for(auto *a: m_lAmbitus) {
    a->paint(painter, index++);
  }

  //нахождение точки по направлению на заданном расстоянии в рад
  SD::punctumRad g0;
  g0.fi = m_fiCurrent;
  g0.lm = m_lmCurrent;
  

  //////////////////////////////////////////////////////////////////
  //нарисовали сетку
  QPen pen;
  pen.setStyle(Qt::DotLine);
  pen.setColor(QColor(Qt::black));
  painter.setPen(pen);

  for(int lon = m_lmMin/M_PI*180; lon<= m_lmMax/M_PI*180; lon+=1) {
    int ll = int(lon)/5*5;
    if(ll!=int(lon))
      continue;
    QPointF p1, p2;
    p1 = fl2xy_grad(m_fiMin/M_PI*180,ll); 
    p2 = fl2xy_grad(m_fiMax/M_PI*180,ll); 
    painter.drawLine(p1,p2);

    ll = int(lon)/20*20;
    if(ll!=int(lon))
      continue;
    painter.drawText(p1.x(), p1.y()/6-10, QObject::tr("%1").arg(ll));
    painter.drawText(p1.x(), 5*p1.y()/6-10, QObject::tr("%1").arg(ll));
  }

  for(int lat = m_fiMin/M_PI*180; lat<= m_fiMax/M_PI*180; lat+=1) {
    int ll = int(lat)/5*5;
    if(ll!=int(lat))
      continue;
    QPointF p1, p2;
    p1 = fl2xy_grad(ll,m_lmMin/M_PI*180); 
    p2 = fl2xy_grad(ll,m_lmMax/M_PI*180); 
    painter.drawLine(p1,p2);

    ll = int(lat)/10*10;
    if(ll!=int(lat))
      continue;
    painter.drawText(p2.x()/6-10,   p2.y()+10, QObject::tr("%1").arg(ll));
    painter.drawText(5*p2.x()/6-10, p2.y()+10, QObject::tr("%1").arg(ll));
    painter.drawText(p2.x()/2-10, p2.y()+10, QObject::tr("%1").arg(ll));

  }
  
  //////////////////////////////////////////////////////////////////
  //траектория по направлению и против
  QVector<int> vAngle;
  vAngle 
    << m_directRad             
    << 180 + m_directRad 
    << m_directRad + m_strobGrad 
    << 180 + m_directRad + m_strobGrad
    << m_directRad - m_strobGrad 
    << 180 + m_directRad - m_strobGrad
    ;

  QPen penDot(Qt::SolidLine);
  penDot.setWidth(2);
  painter.setPen(penDot);
  for(int angle: vAngle) {

    QVector<QPointF> v;
    double lmPrev;

    for(double i=0; i<180; i+=0.2) {
      SD::punctumRad current = SD::operations::toDirection(
        g0, angle/180.*M_PI,i/180.*M_PI);
      if(i==0)
        lmPrev = current.lm;
      current.lm = SD::operations::continuity(lmPrev,current.lm);
      lmPrev = current.lm;
      v << fl2xy(current.fi,current.lm);
    }
    painter.drawPoints(v.constData(),v.size());
  }

  QVector<QPointF> vSegment;
  std::vector<SD::punctumRad> vSegment1;
  SD::operations::buildSector(
    g0,
    m_Lkm, 
    (m_directRad - m_strobGrad)/180.*M_PI,
    (m_directRad + m_strobGrad)/180.*M_PI,
    0.01,
    Re,
    vSegment1
  );
  vSegment.clear();
  for(int i=0; i<vSegment1.size(); ++i)
    vSegment << fl2xy(vSegment1[i].fi,vSegment1[i].lm);
  {
    QPen pen(Qt::yellow);
    pen.setWidth(3);
    painter.setPen(pen);   
    QBrush brush(Qt::Dense4Pattern);
    painter.setBrush(brush);      
    painter.drawPolygon(vSegment.constData(),vSegment.size());
  }

  /////////////////////////////////////////////////////////
  //пересечение с аркой траектория
  //длина арки
  double dbArc = m_Lkm*1000./Re;
  double lmPrev;
  SD::punctumRad currentArc;
  QVector<QPointF> v;
  for(double d=0; d<dbArc; d+=0.001) {
    currentArc = SD::operations::toDirection(
      g0, m_directRad/180.*M_PI,d);
    v << fl2xy(currentArc.fi,currentArc.lm);
  }
  currentArc = SD::operations::toDirection(
    g0, m_directRad/180.*M_PI,dbArc);
  v << fl2xy(currentArc.fi,currentArc.lm);

  {
  QPen pen(Qt::yellow);
  pen.setWidth(3);
  painter.setPen(pen);     
  painter.setBrush(QBrush(Qt::yellow));      
  painter.drawPolyline(v.constData(),v.size());
  }
  /////////////////////////////////////////////////////////
  //пересечение с аркой траектория
  std::vector<SD::punctumOnAmbitus> arcOutputPoints;  
  for(int i=0; i<m_lAmbitus.size(); ++i) {
    m_lAmbitus[i]->a.intersecare(
      m_fiCurrent, m_lmCurrent,     
      currentArc.fi, currentArc.lm,
      &arcOutputPoints
    );
  }

  //
  /////////////////////////////////////////////////////////
  //точка отсчета основаня
  painter.setPen(QColor(Qt::darkGreen));
  painter.setBrush(QBrush(Qt::darkGreen));      
  QPointF point = fl2xy(m_fiCurrent, m_lmCurrent);
  painter.drawEllipse(point, 4, 4);

  /////////////////////////////////////////////////////////
  //точки пересечения с контуром
  qDebug() << "-----------------------------------------------";
  qDebug() << "-----------------------------------------------";
  painter.setPen(Qt::magenta);     
  painter.setBrush(QBrush(Qt::magenta));      
  for(const auto &v: m_vOutputPoints) {
    QPointF point = fl2xy(v.fi, v.lm);
    painter.drawEllipse(point, 4, 4);
    qDebug() << v.idAmbitus << v.idArcus << v.fi << v.lm;
  }

  /////////////////////////////////////////////////////////
  //точки пересечения с аркой
  if(arcOutputPoints.size()) {
    painter.setPen(Qt::cyan);  
    painter.setBrush(QBrush(Qt::cyan));      
    const SD::punctumOnAmbitus &x = arcOutputPoints.at(0);
    int dr = 4;
    QPointF point = fl2xy(x.fi, x.lm);
    painter.drawEllipse(point, dr, dr);
  }
  if(arcOutputPoints.size()>1) {
    painter.setPen(QColor(Qt::black));
    painter.setBrush(QBrush(Qt::black));      
    const SD::punctumOnAmbitus &x = arcOutputPoints.at(arcOutputPoints.size()-1);
    int dr = 4;
    QPointF point = fl2xy(x.fi, x.lm);
    painter.drawEllipse(point, dr, dr);
  }


  painter.setPen(Qt::white);     
  painter.setBrush(QBrush(Qt::white));      

  if(m_numberPunctumByDirection!=-1) {
    QPointF point = fl2xy(m_vOutputPoints.at(m_numberPunctumByDirection).fi, 
      m_vOutputPoints.at(m_numberPunctumByDirection).lm);
    painter.drawEllipse(point, 2, 2);
    painter.drawText(20,20+20,QString::fromLatin1("Dist=%1 angle=%2 strob=%3")
      .arg(m_vOutputPoints.at(m_numberPunctumByDirection).distRad*Re/1000.)
        .arg(m_directRad).arg(m_strobGrad));
  } else {
    painter.drawText(20,20+20,QString::fromLatin1("angle=%1 strob=%2").arg(m_directRad).arg(m_strobGrad));
  }

  int dd=15, aa=30+20, bb=40+20;
  painter.setPen(Qt::darkGreen);     
  painter.setBrush(QBrush(Qt::darkGreen)); 
  painter.drawEllipse(20,aa, 4, 4);
  painter.drawText(30,bb,QObject::tr("- source point"));

  painter.setPen(Qt::magenta);     
  painter.setBrush(QBrush(Qt::magenta)); 
  painter.drawEllipse(20,aa+dd, 4, 4);
  painter.drawText(30,bb+dd,QObject::tr("- intersect direction with contour"));

  painter.setPen(Qt::white);     
  painter.setBrush(QBrush(Qt::white)); 
  painter.drawEllipse(20,aa+2*dd, 4, 4);
  painter.drawText(30,bb+2*dd,QObject::tr("- first intersection by direction"));

  painter.setPen(Qt::cyan);     
  painter.setBrush(QBrush(Qt::cyan)); 
  painter.drawEllipse(20,aa+3*dd, 4, 4);
  painter.drawText(30,bb+3*dd,QObject::tr("- first intersection by direction with limit distance"));

  painter.setPen(Qt::black);     
  painter.setBrush(QBrush(Qt::black)); 
  painter.drawEllipse(20,aa+4*dd, 4, 4);
  painter.drawText(30,bb+4*dd,QObject::tr("- last intersection by direction with limit distance"));

}


//
QPointF Window::fl2xy(double fi, double lm) {
  QPointF p;
  double xM, yM;
  fl2pl(fi, lm, xM, yM);
  p.rx() = m_kx*xM+m_bx;
  p.ry() = m_ky*yM+m_by;
  return p;
}
QPointF Window::fl2xy_grad(double fi, double lm) {
  return fl2xy(fi/180.*M_PI, lm/180.*M_PI);
}
//
void Window::xy2fl(int xPx, int yPx, double &fi, double &lm) {
  QPointF p;
  double xM, yM;
  xM = (xPx - m_bx)/m_kx;
  yM = (yPx - m_by)/m_ky;
  pl2fl(xM, yM, fi, lm);
}

void Window::resizeEvent(QResizeEvent * event) {
  QSize wSize = size();
  resizeData();
}

//
void Window::resizeData() {
  double xMaxPx = size().width();
  double yMaxPx = size().height();
  double yMinPx = 0;
  double xMinPx = 0;

  double xMaxM, yMaxM, yMinM, xMinM;
  double x[4], y[4];

  fl2pl(m_fiMin,m_lmMin,x[0],y[0]);
  fl2pl(m_fiMax,m_lmMax,x[1],y[1]);
  fl2pl(m_fiMin,m_lmMax,x[2],y[2]);
  fl2pl(m_fiMax,m_lmMin,x[3],y[3]);

  xMaxM = qMax(qMax(x[0],x[1]),qMax(x[2],x[3]));
  yMaxM = qMax(qMax(y[0],y[1]),qMax(y[2],y[3]));
  yMinM = qMin(qMin(y[0],y[1]),qMin(y[2],y[3]));
  xMinM = qMin(qMin(x[0],x[1]),qMin(x[2],x[3]));

  m_kx = (xMaxPx-xMinPx)/(xMaxM-xMinM); 
  m_bx = xMinPx-m_kx*xMinM;
  
  m_ky = (yMinPx-yMaxPx)/(yMaxM-yMinM); 
  m_by = yMinPx-m_ky*yMaxM;

  for(auto *a: m_lAmbitus) {
    a->resizeXY(m_kx,m_bx,m_ky,m_by);
  }
  repaint();
}

void ambitusPaint::resizeXY(double kx, double bx, double ky, double by) {
  vXY.clear();
  const std::vector<SD::punctum> &v = a.puncta();
  for(int i=0; i<v.size(); i++) {
    const SD::punctum &p = v.at(i);
    QPointF xyf;
    double xM, yM;
    fl2pl(p.fi, p.lm, xM, yM);
    xyf.rx() = kx*xM + bx;
    xyf.ry() = ky*yM + by;
    vXY << xyf;
  }
}


//
void ambitusPaint::paint(QPainter &painter, int index) {
  //painter.drawPolyline(vXY.constData(),vXY.size());
  painter.drawPoints(vXY.constData(),vXY.size());
}

//
void ambitusPaint::prepare() {
  a.prepare(false);
}

//
void Window::closeEvent(QCloseEvent *event) {
  QSettings settings("sd.ini",  QSettings::IniFormat);  
  settings.setValue("mainWindow/geometry", saveGeometry());
  settings.setValue("mainWindow/windowState", saveState());
  settings.setValue("data/file", m_stAmbitusFile);
  settings.setValue("data/point-fi", m_fiCurrent);
  settings.setValue("data/point-lm", m_lmCurrent);
  QMainWindow::closeEvent(event);
}

void Window::readSettings() {
  QSettings settings("sd.ini",  QSettings::IniFormat);
  restoreGeometry(settings.value("mainWindow/geometry").toByteArray());
  restoreState(settings.value("mainWindow/windowState").toByteArray());
  m_stAmbitusFile = settings.value("data/file").toString();
  m_fiCurrent = settings.value("data/point-fi", m_fiCurrent).toDouble();
  m_lmCurrent = settings.value("data/point-lm", m_lmCurrent).toDouble();
  loadFile();
  recalc();
}


