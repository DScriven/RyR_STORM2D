OBJECTS_DIR   = obj
MOC_DIR       = moc

QT += widgets opengl
CONFIG += c++17

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS
# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    blinkroi.cpp \
    main.cpp \
    cluster.cpp \
    param.cpp \
    setDPI.cpp \
    stormdensity.cpp \
    tiff6qt.cpp \
    viewclusters.cpp

HEADERS += \
    blinkroi.h \
    blinks.h \
    blinkvertex.h \
    cluster.h \
    param.h \
    paramvals.h \
    setDPI.h \
    stormdensity.h \
    tetramer.h \
    tiff6qt.h \
    vertex_data_traits.h \
    viewclusters.h

FORMS = param.ui \
    blinkroi.ui \
    setDPI.ui

RESOURCES = glsl.qrc

unix {
LIBS += -ltiff -lboost_system -lgmp
}
win32 {

LIBS += "C:/dev/CGAL-5.5.1/auxiliary/gmp/lib/gmp.lib"
LIBS += "C:/dev/boost_1_74_0/lib64-msvc-14.2/libboost_system-vc142-mt-x64-1_74.lib"
LIBS +=  "C:/dev/tiff-4.0.10/libtiff/libtiff.lib"
INCLUDEPATH += "C:/dev/tiff-4.0.10/libtiff"
INCLUDEPATH += "C:/dev/CGAL-5.5.1/include" "C:/dev/boost_1_74_0"
}
TARGET = RyR_STORM2D
