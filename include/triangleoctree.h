#ifndef TRIANGLEOCTREE_H
#define TRIANGLEOCTREE_H

#include <octree>
#include <set>
#include <iostream>
#include "triangle.h"
using std::set;

//! The octree node type, values are sets of pointers to Triangles.
typedef octree_node<set<Triangle*> > TSetNode;

//! Octree class for Triangles. The nodes contain a std::set of Triangle pointers, to efficiently search for Triangles.
class TriangleOctree : public octree<set<Triangle*>, 3, std::allocator<set<Triangle* > > >
{
public:
    //! Constructs an empty octree at center with width.
    TriangleOctree( const double* octreeCenter, double octreeWidth );

    //! Adds the triangle to the the smalles possible node in the octree that contains all vertices of the triangle.
    void addTriangle(Triangle *t);

    //! Returns the smallest possible node below the root node that contains all vertices of the given elements.
    //! addChildren determines whether to add children or to stop at leaf nodes.
    //@{
    TSetNode& getNodeForTriangle(const Triangle *t, bool addChildren = true);
    TSetNode& getNodeForTriangle(const Triangle *t, TSetNode &n, const Point &center, const double width, bool addChildren = true);
    TSetNode& getNodeForSphere(const Point &sphereCenter, const double &sphereRadius, bool addChildren = true);
    TSetNode& getNodeForSphere(const Point &sphereCenter, const double &sphereRadius, TSetNode &n, const Point &octantCenter, const double octantWidth, bool addChildren = true);
    TSetNode& getNodeForPointList(List& l, bool addChildren = true );
    TSetNode& getNodeForPointList(const List& l, TSetNode &n, const Point &octantCenter, const double octantWidth, bool addChildren = true );
    //@}
    //! Appends all triangles from and below the node n to the list l.
    void appendTrianglesFromAndBelowNodeToList(TSetNode &n, List &l) const;
    //! Appends all triangles from and above the node n to the list l.
    void appendTrianglesFromAndAboveNodeToList(TSetNode &n, List &l) const;
    //! Appends all triangles from, below and above the node n to the list l.
    void appendTrianglesFromHierarchyToList(TSetNode &n, List &l) const;
    //! Removes all triangles with !isLinked() from and below the node n. Starts from root(), if n not set.
    void removeUnlinkedTriangles( TSetNode *n = NULL );
};
#endif // TRIANGLEOCTREE_H
