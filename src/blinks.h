#ifndef BLINKS_H
#define BLINKS_H
#include <QVector>
#include <QHash>

struct InputBlinks
{
    double x;
    double y;
    double z;
    double dX;
    double dY;
    double dZ;
    double Amplitude;
    int frame;
    int pic_no;
    double R;
};

struct Blinks
{
    double x;
    double y;
    double z;
    double dX;
    double dY;
    double dZ;
    double nnd;
    double Amplitude;
    int xpos;
    int ypos;
    int zpos;
    int pic_no;
};

struct DPoint
{
    double x;
    double y;
};

struct PixelBlinks
{
    int offset;
    int xpos;
    int ypos;
    int zpos;
    int NoBlinks;
};

#endif // BLINKS_H
