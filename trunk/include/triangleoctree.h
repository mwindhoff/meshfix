#ifndef TRIANGLEOCTREE_H
#define TRIANGLEOCTREE_H

#include <octree>
#include <set>
#include <iostream>
#include "triangle.h"
using std::set;
using std::cout;
using std::endl;

//! The octant node type, values are sets of pointers to Triangles.
typedef octree_node<set<Triangle*> > TSetNode;

class TriangleOctree : public octree<set<Triangle*>, 3, std::allocator<set<Triangle* > > >
{
public:
    //! Constructs an empty octree at center with width.
    TriangleOctree( const double* octreeCenter, double octreeWidth );
    //! Adds the triangle to the smallest possible TSetNode, that contains all three vertices of the triangle.
    //! If necessary child nodes are appended.
    void addTriangle(Triangle *t);
    TSetNode& getNodeForTriangle(const Triangle *, TSetNode &n, const Point &center, const double width);
    TSetNode& getNodeForTriangle(const Triangle *);
    TSetNode& getNodeForSphere(const Point &sphereCenter, const double &sphereRadius,
                               TSetNode &n, const Point &octantCenter, const double octantWidth);
    TSetNode& getNodeForSphere(const Point &sphereCenter, const double &sphereRadius);
};
#endif // TRIANGLEOCTREE_H
