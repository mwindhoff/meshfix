#include "triangleoctree.h"
TriangleOctree::TriangleOctree( const double* octreeCenter, double octreeWidth ) :
        octree<List, 3, std::allocator<List> >::octree( octreeCenter, octreeWidth ) {}
TriangleOctree::~TriangleOctree() {
    for(iterator it = this->begin(); it != this->end(); ++it)
        delete(it->valuePtr());
}

void TriangleOctree::addTriangle(Triangle *t) {
    cursor::path p = getPathForTriangle(t);
    p->valuePtr()->appendTail(t);
}
void TriangleOctree::writeTrianglesToLeafs() {
    cursor cs(this);
    Point octantCenter(this->_M_center[0], this->_M_center[1], this->_M_center[2]);
    double octantWidth = this->_M_width;
    writeTrianglesToLeafs(cs.toConstPath(), octantCenter, octantWidth);
}

Point TriangleOctree::path2center(TriangleOctree::cursor::path p, const double *center, double w) {
    Point c(center[0], center[1], center[2]);
    for (int j = 0; j < p._M_indices.size(); j++) {
        int i = p._M_indices[j];
        double w4 = w/4;
        c.x += (i & 1) ? w4 : -1*w4;
        c.y += (i & 2) ? w4 : -1*w4;
        c.z += (i & 4) ? w4 : -1*w4;
        w = w/2;
    }
    return c;
}

void TriangleOctree::writeTrianglesToLeafs(const TriangleOctree::cursor::const_path p, Point octantCenter, double octantWidth) {
    cursor cs(p);
    double w2 = octantWidth/2, w4 = octantWidth/4;
    Point c=octantCenter-Point(w4,w4,w4);
    Triangle *t;
    Node *n;
    if( !cs->is_leaf_node() ) {
        List *l = cs->valuePtr();
        cs.down(0);
        FOREACHVTTRIANGLE(l, t, n) if(triangleIntersectsOctant(t, c, w2))
            cs->valuePtr()->appendHead(t);
        writeTrianglesToLeafs(cs.toConstPath(), c, w2);
        cs.axis_partner(0);
        c.x+=w2;
        FOREACHVTTRIANGLE(l, t, n) if(triangleIntersectsOctant(t, c, w2))
            cs->valuePtr()->appendHead(t);
        writeTrianglesToLeafs(cs.toConstPath(), c, w2);
        cs.axis_partner(1);
        c.y+=w2;
        FOREACHVTTRIANGLE(l, t, n) if(triangleIntersectsOctant(t, c, w2))
            cs->valuePtr()->appendHead(t);
        writeTrianglesToLeafs(cs.toConstPath(), c, w2);
        cs.axis_partner(0);
        c.x-=w2;
        FOREACHVTTRIANGLE(l, t, n) if(triangleIntersectsOctant(t, c, w2))
            cs->valuePtr()->appendHead(t);
        writeTrianglesToLeafs(cs.toConstPath(), c, w2);
        cs.axis_partner(2);
        c.z+=w2;
        FOREACHVTTRIANGLE(l, t, n) if(triangleIntersectsOctant(t, c, w2))
            cs->valuePtr()->appendHead(t);
        writeTrianglesToLeafs(cs.toConstPath(), c, w2);
        cs.axis_partner(0);
        c.x+=w2;
        FOREACHVTTRIANGLE(l, t, n) if(triangleIntersectsOctant(t, c, w2))
            cs->valuePtr()->appendHead(t);
        writeTrianglesToLeafs(cs.toConstPath(), c, w2);
        cs.axis_partner(1);
        c.y-=w2;
        FOREACHVTTRIANGLE(l, t, n) if(triangleIntersectsOctant(t, c, w2))
            cs->valuePtr()->appendHead(t);
        writeTrianglesToLeafs(cs.toConstPath(), c, w2);
        cs.axis_partner(0);
        c.x-=w2;
        FOREACHVTTRIANGLE(l, t, n) if(triangleIntersectsOctant(t, c, w2))
            cs->valuePtr()->appendHead(t);
        writeTrianglesToLeafs(cs.toConstPath(), c, w2);
        // remove the triangles from this node, to save space and because only leafs are of interest
        l->removeNodes();
    }
}

bool TriangleOctree::triangleIntersectsOctant(Triangle *t, Point &octantCenter, double octantWidth) {
    // test wheter the triangle has a point inside the octant
    Point d = (*t->v1())-octantCenter;
    if(fabs(d.x) < octantWidth || fabs(d.y) < octantWidth || fabs(d.z) < octantWidth) return true;
    d = (*t->v2())-octantCenter;
    if(fabs(d.x) < octantWidth || fabs(d.y) < octantWidth || fabs(d.z) < octantWidth) return true;
    d = (*t->v3())-octantCenter;
    if(fabs(d.x) < octantWidth || fabs(d.y) < octantWidth || fabs(d.z) < octantWidth) return true;
    // get the corner points of the octant
    double w2=octantWidth/2;
    Point p1 = octantCenter+Point(-w2,-w2,-w2), // -x -y -z
    p2 = octantCenter+Point(w2,-w2,-w2),        // +x -y -z
    p3 = octantCenter+Point(w2,w2,-w2),         // +x +y -z
    p4 = octantCenter+Point(-w2,w2,-w2),        // -x +y -z
    p5 = octantCenter+Point(-w2,-w2,w2),        // -x -y +z
    p6 = octantCenter+Point(w2,-w2,w2),         // +x -y +z
    p7 = octantCenter+Point(w2,w2,w2),          // +x +y +z
    p8 = octantCenter+Point(-w2,w2,w2);         // -x +y +z
    // test whether the triangle is not completely on the outside of any octant face plane
    Point p11 = (*t->v1())-p1, p12 = (*t->v2())-p1, p13 = (*t->v3())-p1;
    Point p71 = (*t->v1())-p7, p72 = (*t->v2())-p7, p73 = (*t->v3())-p7;
    if(p11.x<0 && p12.x<0 && p13.x<0) return false;
    if(p11.y<0 && p12.y<0 && p13.y<0) return false;
    if(p11.z<0 && p12.z<0 && p13.z<0) return false;
    if(p71.x>0 && p72.x>0 && p73.x>0) return false;
    if(p71.y>0 && p72.y>0 && p73.y>0) return false;
    if(p71.z>0 && p72.z>0 && p73.z>0) return false;
    // intersect triangle plane with the 12 edges, return true if intersection is inside octant and inside the bounding sphere of the triangle
    Point n = t->getNormal(), c = t->getCircleCenter(), i;
    double x,y,z,r2 = (c-(*t->v1())).squaredLength();
    double k = n*(*t->v1());
    if(n.x != 0 ) { // edges in x direction
        x = -(n.y*p1.y+n.z*p1.z-k)/n.x;
        i.setValue(x,p1.y,p1.z);
        if(fabs(x-octantCenter.x) < w2 && i.squaredDistance(&c) < r2 && t->isInside(&i)) return true;
        x = -(n.y*p2.y+n.z*p2.z-k)/n.x;
        i.setValue(x,p2.y,p2.z);
        if(fabs(x-octantCenter.x) < w2 && i.squaredDistance(&c) < r2 && t->isInside(&i)) return true;
        x = -(n.y*p5.y+n.z*p5.z-k)/n.x;
        i.setValue(x,p5.y,p5.z);
        if(fabs(x-octantCenter.x) < w2 && i.squaredDistance(&c) < r2 && t->isInside(&i)) return true;
        x = -(n.y*p6.y+n.z*p6.z-k)/n.x;
        i.setValue(x,p6.y,p6.z);
        if(fabs(x-octantCenter.x) < w2 && i.squaredDistance(&c) < r2 && t->isInside(&i)) return true;
    }
    if(n.y != 0 ) { // edges in y direction
        y = -(n.x*p1.x+n.z*p1.z-k)/n.y;
        i.setValue(y,p1.x,p1.z);
        if(fabs(y-octantCenter.y) < w2 && i.squaredDistance(&c) < r2 && t->isInside(&i)) return true;
        y = -(n.x*p4.x+n.z*p4.z-k)/n.y;
        i.setValue(y,p4.x,p4.z);
        if(fabs(y-octantCenter.y) < w2 && i.squaredDistance(&c) < r2 && t->isInside(&i)) return true;
        y = -(n.x*p5.x+n.z*p5.z-k)/n.y;
        i.setValue(y,p5.x,p5.z);
        if(fabs(y-octantCenter.y) < w2 && i.squaredDistance(&c) < r2 && t->isInside(&i)) return true;
        y = -(n.x*p8.x+n.z*p8.z-k)/n.y;
        i.setValue(y,p8.x,p8.z);
        if(fabs(y-octantCenter.y) < w2 && i.squaredDistance(&c) < r2 && t->isInside(&i)) return true;
    }
    if(n.z != 0 ) { // edges in z direction
        z = -(n.x*p1.x+n.y*p1.y-k)/n.z;
        i.setValue(z,p1.x,p1.y);
        if(fabs(z-octantCenter.z) < w2 && i.squaredDistance(&c) < r2 && t->isInside(&i)) return true;
        z = -(n.x*p2.x+n.y*p2.y-k)/n.z;
        i.setValue(z,p2.x,p2.y);
        if(fabs(z-octantCenter.z) < w2 && i.squaredDistance(&c) < r2 && t->isInside(&i)) return true;
        z = -(n.x*p3.x+n.y*p3.y-k)/n.z;
        i.setValue(z,p3.x,p3.y);
        if(fabs(z-octantCenter.z) < w2 && i.squaredDistance(&c) < r2 && t->isInside(&i)) return true;
        z = -(n.x*p4.x+n.y*p4.y-k)/n.z;
        i.setValue(z,p4.x,p4.y);
        if(fabs(z-octantCenter.z) < w2 && i.squaredDistance(&c) < r2 && t->isInside(&i)) return true;
    }
    return false;
}

TriangleOctree::cursor::path TriangleOctree::getPathForTriangle(const Triangle *t, bool addChildren) {
    cursor cs(this);
    Point octantCenter(this->_M_center[0], this->_M_center[1], this->_M_center[2]);
    double octantWidth = this->_M_width;
    bool b[3];
    while(true) {
        if(cs->is_leaf_node() && !addChildren) break;
        b[0] = t->v1()->x > octantCenter.x;
        b[1] = t->v1()->y > octantCenter.y;
        b[2] = t->v1()->z > octantCenter.z;
        if( b[0] == (t->v2()->x > octantCenter.x) && b[0] == (t->v3()->x > octantCenter.x) &&
            b[1] == (t->v2()->y > octantCenter.y) && b[1] == (t->v3()->y > octantCenter.y) &&
            b[2] == (t->v2()->z > octantCenter.z) && b[2] == (t->v3()->z > octantCenter.z) ) {
            unsigned idx = 0;
            octantCenter -= Point(octantWidth/4, octantWidth/4, octantWidth/4);
            if( b[0] ) {
                idx += 1;
                octantCenter.x += octantWidth/2;
            }
            if( b[1] ) {
                idx += 2;
                octantCenter.y += octantWidth/2;
            }
            if( b[2] ) {
                idx += 4;
                octantCenter.z += octantWidth/2;
            }
            cs->add_children(List()); // if they don't exist yet
            cs.down(idx);
            octantWidth /= 2;
        } else break;
    }
    return cs;
}

TriangleOctree::cursor::path TriangleOctree::getPathForSphere(const Point &sphereCenter, const double &sphereRadius, bool addChildren) {
    cursor cs(this);
    Point octantCenter(this->_M_center[0], this->_M_center[1], this->_M_center[2]);
    double octantWidth = this->_M_width;
    while(true) {
        if(cs->is_leaf_node() && !addChildren) break;
        Point d = (sphereCenter-octantCenter);
        if(MIN(fabs(d.x), MIN(fabs(d.y), fabs(d.z))) > sphereRadius && octantWidth/4 > sphereRadius) {
            unsigned idx = 0;
            octantCenter -= Point(octantWidth/4, octantWidth/4, octantWidth/4);
            if( sphereCenter.x > octantCenter.x ) {
                idx += 1;
                octantCenter.x += octantWidth/2;
            }
            if( sphereCenter.y > octantCenter.y ) {
                idx += 2;
                octantCenter.y += octantWidth/2;
            }
            if( sphereCenter.z > octantCenter.z ) {
                idx += 4;
                octantCenter.z += octantWidth/2;
            }
            cs->add_children(List()); // if they don't exist yet
            cs.down(idx);
            octantWidth /= 2;
        } else break;
    }
    return cs;
}

TriangleOctree::cursor::path TriangleOctree::getPathForPointList(const List& l, bool addChildren) {
    cursor cs(this);
    Point octantCenter(this->_M_center[0], this->_M_center[1], this->_M_center[2]);
    double octantWidth = this->_M_width;
    Point **points = (Point**)l.toArray();
    bool b[3];
    while(cs->is_leaf_node() && !addChildren) {
        b[0] = points[0]->x > octantCenter.x;
        b[1] = points[0]->y > octantCenter.y;
        b[2] = points[0]->z > octantCenter.z;
        unsigned npoints = l.numels();
        for( unsigned i = 1; i < npoints; i++ ) {
            if( b[0] != points[i]->x > octantCenter.x ||
                b[1] != points[i]->y > octantCenter.y ||
                b[2] != points[i]->z > octantCenter.z )
                return cs;
        }
        unsigned idx = 0;
        octantCenter -= Point(octantWidth/4, octantWidth/4, octantWidth/4);
        if( b[0] ) {
            idx += 1;
            octantCenter.x += octantWidth/2;
        }
        if( b[1] ) {
            idx += 2;
            octantCenter.y += octantWidth/2;
        }
        if( b[2] ) {
            idx += 4;
            octantCenter.z += octantWidth/2;
        }
        cs->add_children(List()); // if they don't exist yet
        cs.down(idx);
        octantWidth /= 2;
    }
    return cs;
}

void TriangleOctree::fillTriangleListFromPath(List *l, const TriangleOctree::cursor::const_path p) {
    TriangleOctree::cursor cs(p);
    if(cs->is_leaf_node()) { // append triangles from leaf nodes
        l->appendList(cs->valuePtr());
    } else {
        cs.down(0); // 0 0 0
        fillTriangleListFromPath(l, cs.toConstPath());
        cs.axis_partner(0); // 1 0 0
        fillTriangleListFromPath(l, cs.toConstPath());
        cs.axis_partner(1); // 1 1 0
        fillTriangleListFromPath(l, cs.toConstPath());
        cs.axis_partner(0); // 0 1 0
        fillTriangleListFromPath(l, cs.toConstPath());
        cs.axis_partner(2); // 0 1 1
        fillTriangleListFromPath(l, cs.toConstPath());
        cs.axis_partner(0); // 1 1 1
        fillTriangleListFromPath(l, cs.toConstPath());
        cs.axis_partner(1); // 1 0 1
        fillTriangleListFromPath(l, cs.toConstPath());
        cs.axis_partner(0); // 0 0 1
        fillTriangleListFromPath(l, cs.toConstPath());
    }
}

void TriangleOctree::removeUnlinkedTriangles() {
    iterator it = this->begin(false);
    for( ; it != this->end(false); ++it ) {
        List *l = it->valuePtr();
        List *tmp = new List();
        Triangle *t;
        while( t = (Triangle*) it->valuePtr()->popHead() ) {
            if( t && t->isLinked() )
                tmp->appendTail(t);
        }
        l->appendList(tmp);
    }
}

bool TriangleOctree::isTriangleMaskInPathNode(const TriangleOctree::cursor::const_path &p, char mask) {
    Triangle *t;
    Node *n;
    // check triangles in the node
    FOREACHVTTRIANGLE(p->valuePtr(), t, n) {
        if( t->mask & mask > 0 ) {
            return true;
        }
    }
    return false;
}

bool TriangleOctree::isTriangleMaskInPathDown(TriangleOctree::cursor::path &p, char mask) {
    if(isTriangleMaskInPathNode(p.toConstPath(), mask)) return true;
    if(p->is_leaf_node()) return false;
    // check children
    cursor cs(p.toConstPath());
    cs.down(0); // 0 0 0
    if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
    cs.axis_partner(0); // 1 0 0
    if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
    cs.axis_partner(1); // 1 1 0
    if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
    cs.axis_partner(0); // 0 1 0
    if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
    cs.axis_partner(2); // 0 1 1
    if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
    cs.axis_partner(0); // 1 1 1
    if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
    cs.axis_partner(1); // 1 0 1
    if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
    cs.axis_partner(0); // 0 0 1
    if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
    return false;
}

bool TriangleOctree::findClosestPathContainingTriangleMask(TriangleOctree::cursor::path &p, char mask) {
    // check self and children
    if(isTriangleMaskInPathDown(p, mask)) return true;
    // nothing below found, check nodes above
    cursor cs(p.toConstPath());
    while(cs.level()) {
        cs.axis_partner(0); // 1 0 0
        if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
        cs.axis_partner(1); // 1 1 0
        if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
        cs.axis_partner(0); // 0 1 0
        if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
        cs.axis_partner(2); // 0 1 1
        if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
        cs.axis_partner(0); // 1 1 1
        if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
        cs.axis_partner(1); // 1 0 1
        if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
        cs.axis_partner(0); // 0 0 1
        if(isTriangleMaskInPathDown(cs, mask)) { p = cs; return true; }
        cs.up();            // 0 0 0
        // note that we are not recursing down this one (since we are comming from here):
        if(isTriangleMaskInPathNode(cs.toConstPath(), mask)) { p = cs; return true; }
    }
    return false;
}
