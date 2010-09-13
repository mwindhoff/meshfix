#include "exttrimesh.h"
#include <map>

void ExtTriMesh::initializeOctree() {
    Point b1,b2;
    double bbwidth = getBoundingBox(b1,b2);
    b1=(b1+b2)/2;
    const double a[3] = {b1.x,b1.y,b1.z};
    ot = new TriangleOctree( a, bbwidth );
    Triangle *t2, *t = ((Triangle *)T.head()->data);
    Node *n;
    FOREACHTRIANGLE(t, n) {
        t2 =(Triangle*) n->data;
        ot->addTriangle(t2);
    }
}

/**
 * Joins close or overlapping components of triangles.
 * Algorithm sketch:
 * unmark all triangles
 * for all pairs of components (m,n):
 *  mark2 all triangles of n
 *  mark1 all triangles that are not in n and that are in the hierarchy of n
 *  mark3 all triangles of m that are marked1 and that have a point inside of n
 *  mark3 all triangles that are in the hierarchy of m, have a point inside of m and are marked2 (==part of n)
 *  for all triangles on boundary of m
 *   mark3 all triangles of n within the minAllowedDistance
 *  delete marked3 triangles
 * get the list of boundary loops of m and n
 *  find the closest boundary triangle of n and connect both triangles
 *  skip about 5 triangles
 * fill holes.
 */
int ExtTriMesh::joinCloseOrOverlappingComponents( double minAllowedDistance ) {
    if ( ot == NULL ) initializeOctree();
    Node *n,*m,*nn;
    Triangle *t, *t2;
    unsigned i = 0, j = 0;
    // fill components list
    bool need_unlink = false;
    do {
        need_unlink = false;
        List components;
        fillComponentsList(components);
        // handle all possible pairs of components
        FOREACHNODE(components, n) {
            i++;
            j=0;
            m=n->next();
            if (n->data == m->data && j++ < i) continue;
            List pointList, *hierarchyTriangles, *nl = (List*)n->data, *ml = (List*)m->data;
            // mark3 all triangles of component n, that contain a point inside m
            FOREACHVTTRIANGLE(nl, t, nn) {
                pointList.appendTail(t->v1());
                pointList.appendTail(t->v2());
                pointList.appendTail(t->v3());
                MARK_VISIT2(t);
            }
            ot->appendTrianglesFromHierarchyToList(ot->getNodeForPointList(pointList, false), *hierarchyTriangles);
            bool found = false;
            FOREACHVTTRIANGLE(hierarchyTriangles, t, nn) {
                if( !IS_VISITED2(t) ) { // t not part of component n
                    MARK_VISIT(t);
                    found=true;
                }
            }
            if(found) { // other triangles in the hierarchy of n found, overlap possible
                found = false;
                FOREACHVTTRIANGLE(ml, t, nn) {
                    if( !IS_VISITED(t) && // t was in hierarchy of component n (and could possible overlap with it)
                        ( isPointInComponent((Point*)t->v1(),(List*)n->data) ||
                          isPointInComponent((Point*)t->v2(),(List*)n->data) ||
                          isPointInComponent((Point*)t->v3(),(List*)n->data) )) {
                        MARK_BIT(t,3);
                        need_unlink = true;
                        found = true;
                    }
                }
                if(found) { // triangles of m inside of n found, overlap possible
                    pointList.removeNodes();
                    // mark3 all triangles of component n, that contain a point inside m
                    FOREACHVTTRIANGLE(ml, t, nn) {
                        pointList.appendTail(t->v1());
                        pointList.appendTail(t->v2());
                        pointList.appendTail(t->v3());
                    }
                    FOREACHVTTRIANGLE(hierarchyTriangles, t, nn) UNMARK_VISIT(t);
                    hierarchyTriangles->removeNodes();
                    ot->appendTrianglesFromHierarchyToList(ot->getNodeForPointList(pointList, false), *hierarchyTriangles);
                    FOREACHVTTRIANGLE(hierarchyTriangles, t, nn) {
                        if ( IS_VISITED2(t) && // t was part of component n
                             ( isPointInComponent((Point*)t->v1(),(List*)n->data) ||
                               isPointInComponent((Point*)t->v1(),(List*)n->data) ||
                               isPointInComponent((Point*)t->v1(),(List*)n->data) )) {
                            MARK_BIT(t,3);
                            need_unlink = true;
                        }
                    }
                }
            }
            FOREACHVTTRIANGLE(hierarchyTriangles, t, nn) UNMARK_VISIT(t);
            // mark triangles of n within the minAllowedDistance of any triangle of m for deletion
            FOREACHVTTRIANGLE(ml, t, nn) {
                hierarchyTriangles->removeNodes();
                ot->appendTrianglesFromHierarchyToList(
                        ot->getNodeForSphere(t->getCenter(), minAllowedDistance), *hierarchyTriangles);
                Point circleCenter = t->getCircleCenter();
                FOREACHVTTRIANGLE(hierarchyTriangles, t2, nn) {
                    if (IS_VISITED2(t2) && // is part of component n
                        circleCenter.distance(t2->getCircleCenter()) < minAllowedDistance) {
                        // mark both triangles for deletion
                        MARK_BIT(t2,3);
                        MARK_BIT(t,3);
                        need_unlink = true;
                    }
                }
            }
            // delete marked3 triangles, now both components should have at least one boundary loop
            if (need_unlink) {
                d_boundaries = d_handles = d_shells = 1;
                int deletion_counter = 0;
                nn = nl->head();
                do { // remove from component n
                    t = (Triangle*) nn->data;
                    Node *next = n->prev();
                    if (IS_BIT(t,3)) {
                        t->unlinkEdges();
                        nl->removeCell(n);
                    }
                    n = next;
                } while( n != NULL );
                nn = ml->head();
                do { // remove from component m
                    t = (Triangle*) nn->data;
                    Node *next = n->prev();
                    if (IS_BIT(t,3)) {
                        t->unlinkEdges();
                        ml->removeCell(n);
                    }
                    n = next;
                } while( n != NULL );
                ot->removeUnlinkedTriangles();
                d_boundaries = d_handles = d_shells = 1;
                removeUnlinkedElements();
                // check if both components have a boundary loops and join them
                // TODO: factor 1.5 needed?
                if( joinComponentsCloseBoundaries(nl, ml, minAllowedDistance*1.5 ) ) {
                    // break out of the components loop, reinitialize the components list,
                    // since the number of the components has changed now
                    break;
                }
            }
        }
        // delete components list
        FOREACHNODE(components, n) delete((List *)n->data);
    } while (need_unlink); // until there were no more triangles to delete
    this->unmarkEverything();
    return 0;
}

bool ExtTriMesh::joinComponentsCloseBoundaries(List *nl, List *ml, double maxDistanceToJoin) {
    Vertex *v,*w;
    Node *n, *m, *nn;
    List vertexList, nBoundaryLoops, mBoundaryLoops, *one_loop;
    List **nBoundaryLoopsArray;
    // create list of list boundary loops of component n (= list of vertices)
    getComponentsVertices(*nl, vertexList);
    FOREACHVVVERTEX(((List*)&vertexList), v, n) {
        // find next vertex of an unmarked boundary
        if (IS_VISITED2(v) && !v->isOnBoundary()) continue;
        w = v;
        one_loop = new List;
        do { // mark all vertices at this boundary
            one_loop->appendHead(w);
            MARK_VISIT2(w);
            w = w->nextOnBoundary();
        } while (w != v);
        nBoundaryLoops.appendHead(one_loop);
    }
    FOREACHVVVERTEX(((List*)&vertexList), v, n) UNMARK_VISIT2(v);
    vertexList.removeNodes();
    // create list of list boundary loops of component m (= list of vertices)
    getComponentsVertices(*ml, vertexList);
    FOREACHVVVERTEX(((List*)&vertexList), v, n) {
        // find next vertex of an unmarked boundary
        if (IS_VISITED2(v) && !v->isOnBoundary()) continue;
        w = v;
        one_loop = new List;
        do { // mark all vertices at this boundary
            one_loop->appendHead(w);
            MARK_VISIT2(w);
            w = w->nextOnBoundary();
        } while (w != v);
        mBoundaryLoops.appendHead(one_loop);
    }
    FOREACHVVVERTEX(((List*)&vertexList), v, n) UNMARK_VISIT2(v);

    nBoundaryLoopsArray = (List **)nBoundaryLoops.toArray();
    unsigned nNumLoops = nBoundaryLoops.numels();
    bool allVerticesHaveClosePartner = false;
    double maxDistanceToJoin2 = maxDistanceToJoin*maxDistanceToJoin;
    bool ret = false;
    for (unsigned i = 0; i < nNumLoops; i++) {
        nn = mBoundaryLoops.head();
        do {
            allVerticesHaveClosePartner = false;
            FOREACHVVVERTEX(nBoundaryLoopsArray[i], v, n) {
                bool foundClosePartner = false;
                FOREACHVVVERTEX(((List*)nn->data), w, m) {
                    if( v->squaredDistance(w) < maxDistanceToJoin2 ) {
                        foundClosePartner = true;
                        break;
                    }
                }
                if( !foundClosePartner ) { // no close partner, dont join this pair of boundaries
                    allVerticesHaveClosePartner = false;
                    break;
                }
            }
            if ( allVerticesHaveClosePartner ) {
                joinBoundaryLoops(v, w, true, true, true);
                ret = true;
                mBoundaryLoops.freeCell(nn);
                break;
            }
            nn = nn->next();
        } while ( nn != NULL );
    }
    delete(nBoundaryLoopsArray);
    while ((one_loop=(List *)nBoundaryLoops.popHead())!=NULL) delete one_loop;
    while ((one_loop=(List *)mBoundaryLoops.popHead())!=NULL) delete one_loop;
    return ret;
}

bool ExtTriMesh::isPointInComponent(Point *p, List *component) {
    Node *n;
    Triangle *t, *closestTriangle;
    double mdist = DBL_MAX;
    FOREACHVTTRIANGLE(component, t, n) {
        double dist = MIN(MIN(t->v1()->distance(p), t->v2()->distance(p)), t->v3()->distance(p));
        if( dist < mdist ) {
            mdist = dist;
            closestTriangle = t;
        }
    }
    return closestTriangle->getNormal().getAngle(closestTriangle->getCenter()-p) < M_PI_2;
}

void ExtTriMesh::getComponentsVertices(List &component, List &vertexList) {
    Triangle *t, *t1, *t2, *t3;
    List todo;
    todo.appendHead(component.head());
    while (todo.numels()) {
        t = (Triangle *)todo.popHead();
        if (!IS_VISITED(t)) {
            if(!IS_VISITED(t->v1())) {
                vertexList.appendTail(t->v1());
                MARK_VISIT(t->v1());
            }
            if(!IS_VISITED(t->v2())) {
                vertexList.appendTail(t->v1());
                MARK_VISIT(t->v2());
            }
            if(!IS_VISITED(t->v3())) {
                vertexList.appendTail(t->v1());
                MARK_VISIT(t->v3());
            }
            t1 = t->t1(); t2 = t->t2(); t3 = t->t3();
            if (t1 != NULL && !IS_VISITED(t1)) todo.appendHead(t1);
            if (t2 != NULL && !IS_VISITED(t2)) todo.appendHead(t2);
            if (t3 != NULL && !IS_VISITED(t3)) todo.appendHead(t3);
            MARK_VISIT(t);
        }
    }
    Node *n;
    Vertex *v;
    FOREACHVTTRIANGLE(((List*)&component), t, n) UNMARK_VISIT(t);
    FOREACHVVVERTEX(((List*)&vertexList), v, n) UNMARK_VISIT(v);
}
