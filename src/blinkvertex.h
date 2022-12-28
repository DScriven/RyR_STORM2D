#ifndef BLINKVERTEX_H
#define BLINKVERTEX_H

#include <CGAL/Cartesian.h>

typedef CGAL::Cartesian<double> R;
typedef R::Point_2 Point_2;

struct blinkVertex
{
protected:
    bool   dummy;         // true if the vertex is dummy test point - used to test padding efficiency
    bool   dummyNeighbor; // true if the vertex has at least one dummy point as neighbor - used to test padding efficiency
    Point_2   pos;    // stores particle's position
    double _weight; // weight of each particle
    double _density;// the value of the density will be computed later using the DTFE density interpolation at each vertex position
    int vno;
    
public:
    blinkVertex()
    { 
       _weight = double(1.);
       _density = double(0.);
       vno = -1;
       dummy=false; 
       dummyNeighbor=false; 
    }
    
    inline void setData(blinkVertex &other)
    {
        _weight = other.weight();
        _density = other.density();
        vno=other.VertexNo();
    }


    // do not modify the following
    inline Point_2& position() { return pos;}
    inline const double& position(int const i) { if (i==0)
                                               return pos.x();
                                            else
                                               return pos.y();}
    inline void setPosition(Point_2 p) { pos = p;}
    inline void setVertexNo(int i) { vno = i;}
    inline void setDummy() { dummy=true; dummyNeighbor=true; setDensity(0.); }
    inline void setDummyNeighbor() { dummyNeighbor=true; }
    inline bool isDummy() { return dummy; }
    inline bool hasDummyNeighbor() { return dummyNeighbor; }

    protected:
    
    public:

    // functions to access the weight
    inline int& VertexNo() {return vno;}
    inline double& weight() { return _weight;}
    inline void setWeight(double const w) { _weight = w;}
    // functions to access the density
    inline double& density() { return _density;}
    inline void setDensity(double const d) { _density = d;}
};
#endif
