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

#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <QMainWindow>
#include "../sdist/sdist.h"

struct pointObject {
  QPoint point;
  QRect formular;
  QString text;
  //
  void formularDefault(const QFontMetrics &fm) {
    int dr = 2;
    int pixelsWide = fm.width(text);
    int pixelsHigh = fm.height();
    formular = QRect( point.x() + dr, point.y(), pixelsWide + 2, pixelsHigh);
  }
};

struct ambitusPaint {
  ambitusPaint(int index_) 
    : index(index_) 
    , a(index_)
  {};
  ambitusPaint() 
    : index(0) 
    , a(0)
  {};

  int index;
  SD::ambitus     a;
  QVector<QPointF> vXY; 
  double fiMin, fiMax;
  double lmMin, lmMax;
  //
  void resizeXY(double kx, double bx, double ky, double by);
  //
  void paint(QPainter &painter, int index);
  //
  void prepare();
};

//! [0]
class Window : public QMainWindow
{
    Q_OBJECT

public:
    Window();
    virtual ~Window();
    void readSettings();

private:

  void mouseDoubleClickEvent(QMouseEvent * event);
  void mousePressEvent(QMouseEvent * event);
  void mouseMoveEvent(QMouseEvent * event);
  void paintEvent(QPaintEvent * event);
  void keyReleaseEvent(QKeyEvent * event);
  void keyPressEvent(QKeyEvent * event);
  void resizeEvent(QResizeEvent * event);
  void closeEvent(QCloseEvent *event);

  //  
  void loadFile();
  //  
  void resizeData();
  //
  QPointF fl2xy(double fi, double lm);
  QPointF fl2xy_grad(double fi, double lm);
  void xy2fl(int x, int y, double &fi, double &lm);
  void recalc();

private:
  QList<pointObject> m_lObjects;
  bool m_bDeleteKeyPressed;

  QString m_stAmbitusFile;
  QList<ambitusPaint*> m_lAmbitus;
  std::vector<SD::punctumOnAmbitus> m_vOutputPoints;  

  double m_fiCurrent, m_lmCurrent;
  double m_directRad;
  int m_numberPunctumByDirection;

  double m_fiMin, m_fiMax;
  double m_lmMin, m_lmMax;
  double m_kx, m_bx;
  double m_ky, m_by;
  int m_Lkm;
  int m_Skm;
  int m_strobGrad;
  QVector<SD::punctumRad> m_vMinDistatncePoint;

public slots:
  //  
  void loadDialogFile();

};
//! [0]

#endif // WINDOW_H
