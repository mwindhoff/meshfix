#include "triangleoctree.h"

TriangleOctree::TriangleOctree( const double* octreeCenter, double octreeWidth ) :
        octree<set<Triangle*>, 3, std::allocator<set<Triangle*> > >::octree( octreeCenter, octreeWidth ) {}

void TriangleOctree::addTriangle(Triangle *t) {
    getNodeForTriangle(t).value().insert(t);
}

TSetNode& TriangleOctree::getNodeForTriangle(const Triangle *t, TSetNode &n, const Point &octantCenter, const double octantWidth ) {
    bool b[3];
    b[0] = t->v1()->x > octantCenter.x;
    b[1] = t->v1()->y > octantCenter.y;
    b[2] = t->v1()->z > octantCenter.z;
    if( b[0] == t->v2()->x > octantCenter.x && b[0] == t->v3()->x > octantCenter.x &&
        b[1] == t->v2()->y > octantCenter.y && b[1] == t->v3()->y > octantCenter.y &&
        b[2] == t->v2()->z > octantCenter.z && b[2] == t->v3()->z > octantCenter.z ) {
        unsigned idx = 0;
        Point c = octantCenter-Point(octantWidth/4, octantWidth/4, octantWidth/4);
        if( b[0] ) {
            idx += 1;
            c.x += octantWidth/2;
        }
        if( b[1] ) {
            idx += 2;
            c.y += octantWidth/2;
        }
        if( b[2] ) {
            idx += 4;
            c.z += octantWidth/2;
        }
        n.add_children(); // if they don't exist yet
        return getNodeForTriangle(t, n[idx], c, octantWidth/2);
    }
    return n;
}

TSetNode& TriangleOctree::getNodeForTriangle(const Triangle *t) {
    Point octreeCenter(this->_M_center[0], this->_M_center[1], this->_M_center[2]);
    return getNodeForTriangle(t, *this->root(), octreeCenter, this->_M_size);
}

TSetNode& TriangleOctree::getNodeForSphere(const Point &sphereCenter, const double &sphereRadius,
                                           TSetNode &n, const Point &octantCenter, const double octantWidth) {
    if( n.num_children() > 0 && sphereCenter.distance(octantCenter) > sphereRadius ) {
        unsigned idx = 0;
        Point c = octantCenter-Point(octantWidth/4, octantWidth/4, octantWidth/4);
        if( sphereCenter.x > octantCenter.x ) {
            idx += 1;
            c.x += octantWidth/2;
        }
        if( sphereCenter.y > octantCenter.y ) {
            idx += 2;
            c.y += octantWidth/2;
        }
        if( sphereCenter.z > octantCenter.z ) {
            idx += 4;
            c.z += octantWidth/2;
        }
        return getNodeForSphere(sphereCenter, sphereRadius, n[idx], c, octantWidth/2);
    }
    return n;
}

TSetNode& TriangleOctree::getNodeForSphere(const Point &sphereCenter, const double &sphereRadius) {
    Point octreeCenter(this->_M_center[0], this->_M_center[1], this->_M_center[2]);
    return getNodeForSphere(sphereCenter, sphereRadius, *this->root(), octreeCenter, this->_M_size);
}
