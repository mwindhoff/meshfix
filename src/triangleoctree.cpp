#include "triangleoctree.h"

TriangleOctree::TriangleOctree( const double* center, double size ) :
        octree<set<Triangle*>, 3, std::allocator<set<Triangle*> > >::octree( center, size ) {}

void TriangleOctree::addTriangle(Triangle *t) {
    getNodeForTriangle(t).value().insert(t);
}

octree_node<set<Triangle*> >& TriangleOctree::getNodeForTriangle(Triangle *t, octree_node<set<Triangle*> >& n, const double *center, double width ) {
    bool b[3][3];
    b[0][0]= t->v1()->x > center[0];
    b[0][1]= t->v1()->y > center[1];
    b[0][2]= t->v1()->z > center[2];
    b[1][0]= t->v2()->x > center[0];
    b[1][1]= t->v2()->y > center[1];
    b[1][2]= t->v2()->z > center[2];
    b[2][0]= t->v3()->x > center[0];
    b[2][1]= t->v3()->y > center[1];
    b[2][2]= t->v3()->z > center[2];
    if(    (b[0][0] == b[1][0] == b[2][0])
        && (b[0][1] == b[1][1] == b[2][1])
        && (b[0][2] == b[1][2] == b[2][2])) {
        unsigned idx = 0;
        double c[3] = {center[0], center[1], center[2]};
        for ( int i = 0; i < 3; i++) {
            if (b[i][0]) {
                idx+=(1<<i);
                c[i] += width/2;
            } else
                c[i] -= width/2;
        }
        n.add_children(); // if they don't exist yet
        return getNodeForTriangle(&*t, n[idx], (const double*) c, width/2);
    }
    return n;
}

octree_node<set<Triangle*> >& TriangleOctree::getNodeForTriangle(Triangle *t) {
    return getNodeForTriangle(&*t, *this->root(), this->center(), this->size());
}
