#include "exttrimesh.h"
#include <map>

void ExtTriMesh::initializeOctree() {
    Point b1,b2;
    double bbwidth = getBoundingBox(b1,b2);
    b1=(b1+b2)/2;
    const double a[3] = {b1.x,b1.y,b1.z};
    ot = new TriangleOctree( a, bbwidth );
    Triangle *t;
    Node *n;
    FOREACHTRIANGLE(t, n) {
        ot->addTriangle(t);
    }
//    printf("nT: %d\n", T.numels());
//    TriangleOctree::iterator it = ot->begin(false);
//    int sum = 0, counter = 0;
//    for ( ; it != ot->end(false); ++it ) {
//        sum += it->value().numels();
//        if( it->value().numels() > 0 ) {
//            FOREACHVTTRIANGLE(it->valuePtr(), t, n) {
//                if( t->getCircleCenter().distance(t->v1()) > bbwidth/(1<<(it->getLevel())) )
//                    counter++;
//            }
//        }
//    }
//    printf( "sum: %d, counter: %d, badCounter: %d\n", sum, counter);
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
    JMesh::info("joinCloseOrOverlappingComponents():\n");
    JMesh::info(" Initializing octree\n");
    initializeOctree();
    Node *n, *m, *nn, *mm;
    Triangle *t, *t2;
    unsigned i = 0, j = 0;
    // fill components list
    bool need_unlink = false;
    do {
        i = 0;
        need_unlink = false;
        List *components = getComponentsList();
        JMesh::info(" Number of components: %d\n", components->numels() );
        // handle all possible pairs of components
        FOREACHNODE(*components, n) {
            m=n->next();
            j=i+1;
            do {
                if (m == NULL || n == NULL || n->data == m->data) break;
                JMesh::info(" Checking components %d,%d for overlap.\n", i, j);
                List pointList, *nl = (List*)n->data, *ml = (List*)m->data;
                // mark1 triangles of component n
                FOREACHVTTRIANGLE(nl, t, nn) {
                    pointList.appendTail(t->v1());
                    pointList.appendTail(t->v2());
                    pointList.appendTail(t->v3());
                    MARK_BIT(t,1);
                }
                printf("13 %d\n", pointList.numels());
                List *nHierarchyTriangles = ot->getTriangleListFromPath(
                        ot->getPathForPointList(pointList, false).toConstPath());
                bool found = false;
                printf("14 %d\n", nHierarchyTriangles->numels());
                // mark0 triangles in hierarchy of n, that are not part of n
                FOREACHVTTRIANGLE(nHierarchyTriangles, t, nn) {
                    if( !IS_BIT(t,1) ) { // t not part of component n
                        MARK_BIT(t,0);
                        found=true;
                    }
                }
                printf("15\n");
                if(found) { // other triangles in the hierarchy of n found, overlap possible
                    JMesh::info(" Other triangles in the hierarchy of %d found, overlap possible.\n",i);
                    found = false;
                    // mark2 all triangles of component n, that contain a point inside m
                    FOREACHVTTRIANGLE(ml, t, nn) {
                        if( IS_BIT(t,0) && // t was in hierarchy of component n (and could possible overlap with it)
                            ( isPointInComponent(t->v1(),1) ||
                              isPointInComponent(t->v2(),1) ||
                              isPointInComponent(t->v3(),1) )) {
                            printf("10\n");
                            MARK_BIT(t,2);
                            need_unlink = true;
                            found = true;
                        }
                    }
                    printf("12\n");
                    if(found) { // triangles of m inside of n found, reverse overlap possible
                        JMesh::info(" Triangles of %d inside of %d found, searching for triangles of %d in %d.", j, i, i, j);
                        pointList.removeNodes();
                        // mark3 triangles of m
                        FOREACHVTTRIANGLE(ml, t, nn) {
                            pointList.appendTail(t->v1());
                            pointList.appendTail(t->v2());
                            pointList.appendTail(t->v3());
                            MARK_BIT(t,3);
                        }
                        List *mHierarchyTriangles = ot->getTriangleListFromPath(
                                ot->getPathForPointList(pointList, false).toConstPath());
                        // mark2 all triangles of component n, that contain a point inside m
                        FOREACHVTTRIANGLE(mHierarchyTriangles, t, nn) {
                            if ( IS_BIT(t,1) && // t was part of component n
                                 ( isPointInComponent((Point*)t->v1(),3) ||
                                   isPointInComponent((Point*)t->v2(),3) ||
                                   isPointInComponent((Point*)t->v3(),3) )) {
                                MARK_BIT(t,2);
                                need_unlink = true;
                            }
                        }
                        mHierarchyTriangles->removeNodes();
                        // unmark3 triangles of m again
                        FOREACHVTTRIANGLE(ml, t, nn) UNMARK_BIT(t,3);
                    } else JMesh::info(" No triangles of %d in %d found.\n", i, j);
                }
                printf("16\n");
                // unmark0 triangles in hierarchy of n again, that are not part of n
                FOREACHVTTRIANGLE(nHierarchyTriangles, t, nn) UNMARK_BIT(t,0);
                // mark triangles of n within the minAllowedDistance of any triangle of m for deletion
                double minAllowedDistance2 = minAllowedDistance*minAllowedDistance;
                FOREACHVTTRIANGLE(ml, t, nn) {
                    nHierarchyTriangles->removeNodes();
                    nHierarchyTriangles = ot->getTriangleListFromPath(
                            ot->getPathForSphere(t->getCenter(), minAllowedDistance).toConstPath());
                    Point circleCenter = t->getCircleCenter();
                    FOREACHVTTRIANGLE(nHierarchyTriangles, t2, mm) {
                        if (IS_BIT(t2,1) && // is part of component n
                            circleCenter.squaredDistance(t2->getCircleCenter()) < minAllowedDistance2) {
                            // mark both triangles for deletion
                            MARK_BIT(t2,2);
                            MARK_BIT(t,2);
                            need_unlink = true;
                        }
                    }
                }
                nHierarchyTriangles->removeNodes();
                // delete marked2 triangles, now both components should have at least one boundary loop
                if (need_unlink) {
                    d_boundaries = d_handles = d_shells = 1;
                    int deletion_counter = 0;
                    List tmp;
                    while( t = (Triangle*) nl->popHead() ) {
                        if (IS_BIT(t,2)) {
                            t->unlinkEdges();
                            nl->removeCell(nn);
                            deletion_counter++;
                        }
                    }
                    nn = nl->head();
                    do { // remove from component n
                        t = (Triangle*) nn->data;
                        Node *next = nn->prev();
                        if (IS_BIT(t,2)) {
                            t->unlinkEdges();
                            nl->removeCell(nn);
                            deletion_counter++;
                        }
                        nn = next;
                    } while( nn != NULL );
                    nn = ml->head();
                    do { // remove from component m
                        t = (Triangle*) nn->data;
                        Node *next = nn->prev();
                        if (IS_BIT(t,2)) {
                            t->unlinkEdges();
                            ml->removeCell(nn);
                            deletion_counter++;
                        }
                        nn = next;
                    } while( nn != NULL );
                    JMesh::info(" Deleted %d triangles\n", deletion_counter);
                    ot->removeUnlinkedTriangles();
                    d_boundaries = d_handles = d_shells = 1;
                    removeUnlinkedElements();
                    // check if both components have a boundary loops and join them
                    // TODO: factor 1.5 needed?
                    JMesh::info(" Joining close boundaries of %d,%d\n", i,j);
                    if( joinComponentsCloseBoundaries(nl, ml, minAllowedDistance*1.5 ) ) {
                        // break out of the components loop, reinitialize the components list,
                        // since the number of the components has changed now
                        break;
                    }
                }
                m = m->next();
                j++;
            } while( m != NULL );
            i++;
        }
        // delete components list
        FOREACHNODE(*components, n) delete((List *)n->data);
    } while (need_unlink); // until there were no more triangles to delete
    this->unmarkEverything();
    return 0;
}

bool ExtTriMesh::joinComponentsCloseBoundaries(List *nl, List *ml, double maxDistanceToJoin) {
    Vertex *v,*w;
    Node *n, *m, *nn;
    List *vertexList, nBoundaryLoops, mBoundaryLoops, *one_loop;
    List **nBoundaryLoopsArray;
    // create list of list boundary loops of component n (= list of vertices)
    vertexList = getComponentsVertices(nl);
    FOREACHVVVERTEX(vertexList, v, n) {
        // find next vertex of an unmarked boundary
        if (IS_BIT(v,1) && !v->isOnBoundary()) continue;
        w = v;
        one_loop = new List;
        do { // mark all vertices at this boundary
            one_loop->appendHead(w);
            MARK_BIT(w,1);
            w = w->nextOnBoundary();
        } while (w != v);
        nBoundaryLoops.appendHead(one_loop);
    }
    FOREACHVVVERTEX(vertexList, v, n) UNMARK_BIT(v,1);
    vertexList->removeNodes();
    // create list of list boundary loops of component m (= list of vertices)
    vertexList = getComponentsVertices(ml);
    FOREACHVVVERTEX(vertexList, v, n) {
        // find next vertex of an unmarked boundary
        if (IS_BIT(v,1) && !v->isOnBoundary()) continue;
        w = v;
        one_loop = new List;
        do { // mark all vertices at this boundary
            one_loop->appendHead(w);
            MARK_BIT(w,1);
            w = w->nextOnBoundary();
        } while (w != v);
        mBoundaryLoops.appendHead(one_loop);
    }
    FOREACHVVVERTEX(vertexList, v, n) UNMARK_BIT(v,1);

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

bool ExtTriMesh::isPointInComponent(Point *p, char bit) {
    TriangleOctree::cursor::path n = ot->getPathForSphere(*p, 1.0, false);
    if ( ot->findClosestPathContainingTriangleMask(n, 1<<bit) ) {
        List *l = ot->getTriangleListFromPath(n.toConstPath());
        double mdist = DBL_MAX;
        bool found = false;
        Triangle *t, *closestTriangle;
        Node *nn;
        int num = l->numels();
        FOREACHVTTRIANGLE(l, t, nn) {
            if( (t->mask & mask) > 0 ) { // part of component
                double dist = t->getCenter().squaredDistance(p);
                if( dist < mdist ) {
                    mdist = dist;
                    closestTriangle = t;
                    found = true;
                }
            }
        }
        l->removeNodes();
        if(!found)
            return false;
        return closestTriangle->getNormal().getAngle(closestTriangle->getCenter()-p) < M_PI_2;
    }
    return false;
}

List* ExtTriMesh::getComponentsVertices(List *component) {
    Triangle *t, *t1, *t2, *t3;
    List todo;
    List *vertexList = new List();
    todo.appendHead(component->head());
    while (todo.numels()) {
        t = (Triangle *)todo.popHead();
        if (!IS_BIT(t,0)) {
            if(!IS_BIT(t->v1(),0)) {
                vertexList->appendTail(t->v1());
                MARK_BIT(t->v1(),0);
            }
            if(!IS_BIT(t->v2(),0)) {
                vertexList->appendTail(t->v1());
                MARK_BIT(t->v2(),0);
            }
            if(!IS_BIT(t->v3(),0)) {
                vertexList->appendTail(t->v1());
                MARK_BIT(t->v3(),0);
            }
            t1 = t->t1(); t2 = t->t2(); t3 = t->t3();
            if (t1 != NULL && !IS_BIT(t1,0)) todo.appendHead(t1);
            if (t2 != NULL && !IS_BIT(t2,0)) todo.appendHead(t2);
            if (t3 != NULL && !IS_BIT(t3,0)) todo.appendHead(t3);
            MARK_BIT(t,0);
        }
    }
    Node *n;
    Vertex *v;
    FOREACHVTTRIANGLE(component, t, n) UNMARK_BIT(t,0);
    FOREACHVVVERTEX(vertexList, v, n) UNMARK_BIT(v,0);
    return vertexList;
}
