
#include "exttrimesh.h"
#include <map>

int ExtTriMesh::initializeOctree() {
    Triangle *t = ((Triangle *)T.head()->data);
    const double a[3] = {.0,.0,.0};
    ot = new octree( a, 128., t);
    Node *n;
    FOREACHTRIANGLE(t, n) {
        n->data
    }
}

int ExtTriMesh::joinCloseOrOverlappingComponents( float min_allowed_distance ) {
    Node *n,*m;
    std::map<const unsigned, const List*> sizeListMap;
    Triangle *t;
    int nt = 0;
    // fill components list
    List components;
    components = getComponents();
    FOREACHNODE(components, n)
        sizeListMap[((unsigned)((List *)n->data)->numels())]=(List *)n->data;
    FOREACHTRIANGLE(t, n) UNMARK_VISIT2(t);
    // handle all possible pairs of components
    FOREACHNODE(components, n) {
        FOREACHNODE(components, m) {
            if (n == m) break;
            // mark all triangles of component n, that contain a point inside m
            FOREACHTRIANGLE(t, n) {
                ;
            }
        }
    }

    // delete components list
    FOREACHNODE(components, n) delete((List *)n->data);
    // if there are components that were unlinked
    if (nt) {
        d_boundaries = d_handles = d_shells = 1;
        removeUnlinkedElements();
        return 1;
    }
    return 0;
}


#include "exttrimesh.h"
#include <map>

int ExtTriMesh::initializeOctree() {
    Triangle *t = ((Triangle *)T.head()->data);
    const double a[3] = {.0,.0,.0};
    ot = new octree( a, 128., t);
    Node *n;
    FOREACHTRIANGLE(t, n) {
        n->data
    }
}

int ExtTriMesh::joinCloseOrOverlappingComponents( float min_allowed_distance ) {
    Node *n,*m;
    std::map<const unsigned, const List*> sizeListMap;
    Triangle *t;
    int nt = 0;
    // fill components list
    List components;
    components = getComponents();
    FOREACHNODE(components, n)
        sizeListMap[((unsigned)((List *)n->data)->numels())]=(List *)n->data;
    FOREACHTRIANGLE(t, n) UNMARK_VISIT2(t);
    // handle all possible pairs of components
    FOREACHNODE(components, n) {
        FOREACHNODE(components, m) {
            if (n == m) break;
            // mark all triangles of component n, that contain a point inside m
            FOREACHTRIANGLE(t, n) {
                ;
            }
        }
    }

    // delete components list
    FOREACHNODE(components, n) delete((List *)n->data);
    // if there are components that were unlinked
    if (nt) {
        d_boundaries = d_handles = d_shells = 1;
        removeUnlinkedElements();
        return 1;
    }
    return 0;
}

