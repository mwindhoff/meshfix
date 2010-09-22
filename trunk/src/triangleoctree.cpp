#include "triangleoctree.h"
TriangleOctree::TriangleOctree( const double* octreeCenter, double octreeWidth ) :
        octree<List, 3, std::allocator<List> >::octree( octreeCenter, octreeWidth ) {}

void TriangleOctree::addTriangle(Triangle *t) {
    cursor::path p = getPathForTriangle(t);
    p->valuePtr()->appendTail(t);
}

TriangleOctree::cursor::path TriangleOctree::getPathForTriangle(const Triangle *t, bool addChildren) {
    cursor cs(this);
    Point octantCenter(this->_M_center[0], this->_M_center[1], this->_M_center[2]);
    double octantWidth = this->_M_size;
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
    double octantWidth = this->_M_size;
    while(true) {
        if(cs->is_leaf_node() && !addChildren) break;
        if( sphereCenter.distance(octantCenter) > sphereRadius ) {
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

TriangleOctree::cursor::path TriangleOctree::getPathForPointList(List& l, bool addChildren) {
    cursor cs(this);
    Point octantCenter(this->_M_center[0], this->_M_center[1], this->_M_center[2]);
    double octantWidth = this->_M_size;
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

List* TriangleOctree::getTriangleListFromPathDown(const TriangleOctree::cursor::const_path p) {
    TriangleOctree::cursor cs(p);
    List *l = new List();
    if(!cs->is_leaf_node()) {
        cs.down(0); // 0 0 0
        l->joinTailList(getTriangleListFromPathDown(cs.toConstPath()));
        cs.axis_partner(0); // 1 0 0
        l->joinTailList(getTriangleListFromPathDown(cs.toConstPath()));
        cs.axis_partner(1); // 1 1 0
        l->joinTailList(getTriangleListFromPathDown(cs.toConstPath()));
        cs.axis_partner(0); // 0 1 0
        l->joinTailList(getTriangleListFromPathDown(cs.toConstPath()));
        cs.axis_partner(2); // 0 1 1
        l->joinTailList(getTriangleListFromPathDown(cs.toConstPath()));
        cs.axis_partner(0); // 1 1 1
        l->joinTailList(getTriangleListFromPathDown(cs.toConstPath()));
        cs.axis_partner(1); // 1 0 1
        l->joinTailList(getTriangleListFromPathDown(cs.toConstPath()));
        cs.axis_partner(0); // 0 0 1
        l->joinTailList(getTriangleListFromPathDown(cs.toConstPath()));
    }
    return l;
}

List* TriangleOctree::getTriangleListFromPathUp(const TriangleOctree::cursor::const_path p) {
    cursor cs(p);
    List *l = cs->valuePtr();
    while(cs.level()) {
        cs.up();
        l->joinTailList(cs->valuePtr());
    }
    return l;
}

List* TriangleOctree::getTriangleListFromPath(const TriangleOctree::cursor::const_path p) {
    List *l = getTriangleListFromPathUp(p);
    l->joinTailList(getTriangleListFromPathDown(p));
    return l;
}

void TriangleOctree::removeUnlinkedTriangles() {
    iterator it = this->begin(false);
    for( ; it != this->end(false); ++it ) {
        List *l = it->valuePtr();
        List *tmp = new List();
        Triangle *t;
        while( t = (Triangle*) it->value().popHead() ) {
            if( t && t->isLinked() )
                tmp->appendTail(t);
        }
        l->joinTailList(tmp);
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
