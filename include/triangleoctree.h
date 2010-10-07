#ifndef TRIANGLEOCTREE_H
#define TRIANGLEOCTREE_H

#include <octree>
#include <set>
#include "triangle.h"

using std::set;

typedef octree_node<List, 3, std::allocator<List> > TriangleOctreeNode;
//! Octree class for Triangles. The nodes contain a std::set of Triangle pointers, to efficiently search for Triangles.
class TriangleOctree : public octree<List, 3, std::allocator<List> >
{
public:
    //! Constructs an empty octree at center with width.
    TriangleOctree( const double* octreeCenter, double octreeWidth );
    ~TriangleOctree();
    //! Adds the triangle to the the smalles possible node in the octree that contains all vertices of the triangle.
    void addTriangle(Triangle *t);
    //! Appends triangles of non leaf nodes to the children, which they intersect. Afterwards all triangles are in the leaf nodes.
    void writeTrianglesToLeafs();
    void writeTrianglesToLeafs(const cursor::const_path p, Point octantCenter, double octantWidth);
    //! Tests whether the Triangle intersects with the octant.
    bool triangleIntersectsOctant(Triangle *t, Point &octantCenter, double octantWidth);
    //! Returns the smallest possible node below the root node that contains all vertices of the given elements.
    //! addChildren determines whether to add children or to stop at leaf nodes.
    //@{
    cursor::path getPathForTriangle(const Triangle *t, bool addChildren = true);
    cursor::path getPathForSphere(const Point &sphereCenter, const double &sphereRadius, bool addChildren = true);
    cursor::path getPathForPointList(const List &l, bool addChildren = true );
    //@}
    //! Returns a list of all triangles from below the node indicated by the path.
    void fillTriangleListFromPath(List *l, const TriangleOctree::cursor::const_path p);
    //! Removes all triangles with !isLinked() from and below the node n. Starts from root(), if n not set.
    void removeUnlinkedTriangles();
    //! Checks, whether a triangle with (t->mask & mask) > 0 is contained in the node indicated by the path.
    bool isTriangleMaskInPathNode(const TriangleOctree::cursor::const_path &p, char mask);
    //! Checks, whether a triangle with (t->mask & mask) > 0 is contained below the node indicated by the path.
    bool isTriangleMaskInPathDown(TriangleOctree::cursor::path &p, char mask);
    //! Finds the path to the first node of the hierarchy with (t->mask & mask) > 0, which encloses the inital path
    //! and the triangle. Means that if such a triangle is found below the path p, p remains unchanged.
    bool findClosestPathContainingTriangleMask(TriangleOctree::cursor::path &p, char mask);

    //! debug stuff
    Point path2center(TriangleOctree::cursor::path p, const double *center, double w);
};
#endif // TRIANGLEOCTREE_H
