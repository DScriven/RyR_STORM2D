#ifndef BLINKVERTEX_TRAITS_H
#define BLINKVERTEX_TRAITS_H
#include "blinkvertex.h"

struct blinkVertex_less_x
{
    bool operator()(blinkVertex p, blinkVertex q) const
    { return (p.position(0) < q.position(0)); }
};
struct blinkVertex_less_y
{
    bool operator()(blinkVertex p, blinkVertex q) const
    { return (p.position(1) < q.position(1)); }
};
struct blinkVertex_sort_traits
{
    typedef blinkVertex Point_2;
    typedef blinkVertex_less_x Less_x_2;
    typedef blinkVertex_less_y Less_y_2;

    Less_x_2 less_x_2_object() const
    { return Less_x_2(); }
    Less_y_2 less_y_2_object() const
    { return Less_y_2(); }
    blinkVertex_sort_traits() {}
};


// Compares two Vertices to be able to sort them according to positions

bool operator < ( blinkVertex p1, blinkVertex p2)
{
    return ((p1.position(0) < p2.position(0)) && (p1.position(1) < p2.position(1)));
}

/*
template <typename T>
bool sameVertex(T p1, T p2)
{
    return (p1.position(0)==p2.position(0) && p1.position(1)==p2.position(1));
}
*/
#endif
// blinkVertex_TRAITS_H
