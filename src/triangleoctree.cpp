#include "triangleoctree.h"

TriangleOctree::TriangleOctree( const double* octreeCenter, double octreeWidth ) :
        octree<set<Triangle*>, 3, std::allocator<set<Triangle*> > >::octree( octreeCenter, octreeWidth ) {}

void TriangleOctree::addTriangle(Triangle *t) {
    getNodeForTriangle(t).value().insert(t);
}

TSetNode& TriangleOctree::getNodeForTriangle(const Triangle *t, bool addChildren) {
    Point octreeCenter(this->_M_center[0], this->_M_center[1], this->_M_center[2]);
    return getNodeForTriangle(t, *this->root(), octreeCenter, this->_M_size);
}

TSetNode& TriangleOctree::getNodeForTriangle(const Triangle *t, TSetNode &n, const Point &octantCenter, const double octantWidth, bool addChildren ) {
    if( !addChildren && n.is_leaf_node() ) return n;
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

TSetNode& TriangleOctree::getNodeForSphere(const Point &sphereCenter, const double &sphereRadius, bool addChildren ) {
    Point octreeCenter(this->_M_center[0], this->_M_center[1], this->_M_center[2]);
    return getNodeForSphere(sphereCenter, sphereRadius, *this->root(), octreeCenter, this->_M_size, addChildren );
}

TSetNode& TriangleOctree::getNodeForSphere(const Point &sphereCenter, const double &sphereRadius,
                                           TSetNode &n, const Point &octantCenter, const double octantWidth, bool addChildren) {
    if( !addChildren && n.is_leaf_node() ) return n;
    if( sphereCenter.distance(octantCenter) > sphereRadius ) {
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
        n.add_children(); // if they don't exist yet
        return getNodeForSphere(sphereCenter, sphereRadius, n[idx], c, octantWidth/2);
    }
    return n;
}

TSetNode& TriangleOctree::getNodeForPointList(List& l, bool addChildren) {
    Point octreeCenter(this->_M_center[0], this->_M_center[1], this->_M_center[2]);
    return getNodeForPointList(l, *this->root(), octreeCenter, this->_M_size, addChildren );
}

TSetNode& TriangleOctree::getNodeForPointList(const List& l, TSetNode &n, const Point &octantCenter, const double octantWidth, bool addChildren ) {
    if( !addChildren && n.is_leaf_node() ) return n;
    Point **points = (Point**)l.toArray();
    bool b[3];
    b[0] = points[0]->x > octantCenter.x;
    b[1] = points[0]->y > octantCenter.y;
    b[2] = points[0]->z > octantCenter.z;
    unsigned npoints = l.numels();
    for( unsigned i = 1; i < npoints; i++ ) {
        if( b[0] != points[i]->x > octantCenter.x ||
            b[1] != points[i]->y > octantCenter.y ||
            b[2] != points[i]->z > octantCenter.z )
            return n;
    }
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
    return getNodeForPointList(l, n[idx], c, octantWidth/2);
}


void TriangleOctree::appendTrianglesFromAndBelowNodeToList(TSetNode &n, List &l) const {
	if ( n.value().size() > 0 ) {
		set<Triangle*>::iterator it;
		for ( it = n.value().begin(); it != n.value().end(); it++ )
			l.appendTail(*it);
	}
	if ( n.num_children() == 8 )
		for( unsigned i = 0; i < 8; i++ )
			appendTrianglesFromAndBelowNodeToList(n[i], l);
}

void TriangleOctree::appendTrianglesFromAndAboveNodeToList(TSetNode &n, List &l) const {
	if ( n.value().size() > 0 ) {
		set<Triangle*>::iterator it;
		for ( it = n.value().begin(); it != n.value().end(); it++ )
			l.appendTail(*it);
	}
    if ( n._M_parent != NULL && n._M_parent != this->_M_root )
        appendTrianglesFromAndAboveNodeToList(*n._M_parent, l);
}

void TriangleOctree::appendTrianglesFromHierarchyToList(TSetNode &n, List &l) const {
    appendTrianglesFromAndBelowNodeToList(n, l);
    if ( n._M_parent != NULL && n._M_parent != this->_M_root )
        appendTrianglesFromAndAboveNodeToList(*n._M_parent, l);
}

void TriangleOctree::removeUnlinkedTriangles( TSetNode *n ) {
    if ( n == NULL ) n=this->root();
    if ( !n->is_leaf_node() )
        for ( unsigned i = 0; i < 8; i++ )
            removeUnlinkedTriangles(&n[i]);
    set<Triangle*>::iterator it;
    for ( it = n->value().begin(); it != n->value().end(); it++ );
        if( !(*it)->isLinked() )
            n->value().erase(it);
}
