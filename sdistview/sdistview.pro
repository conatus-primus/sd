QT += widgets

LIBS_PATH =  C:/OSGeo4W/lib

HEADERS       = window.h
SOURCES       = main.cpp \
				window.cpp\
				../sdist/sdist.cpp

LIBS += -L$${LIBS_PATH}/
LIBS += -lproj_i




