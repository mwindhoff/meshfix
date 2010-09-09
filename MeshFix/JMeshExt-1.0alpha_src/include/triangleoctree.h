#ifndef TRIANGLEOCTREE_H
#define TRIANGLEOCTREE_H

#include <octree>
#include <set>
#include <iostream>
#include "triangle.h"
using std::set;
using std::cout;
using std::endl;

class TriangleOctree : public octree<set<Triangle*>, 3, std::allocator<set<Triangle* > > >
{
public:
    TriangleOctree( const double* center, double size );
    void addTriangle(Triangle *t);
    octree_node<set<Triangle*> >& getNodeForTriangle(Triangle *, octree_node<set<Triangle*> > &n, const double *center, double width);
    octree_node<set<Triangle*> >& getNodeForTriangle(Triangle *);
};
#endif // TRIANGLEOCTREE_H
