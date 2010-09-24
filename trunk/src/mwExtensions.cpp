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
    struct {
        List *triangles;
        List *hierarchyTriangles;
        Point center;
        double inSphereRadius2;
        double outSphereRadius2;
    } cn, cm; // component structs
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
                List pointList;
                cn.triangles = (List*)n->data;
                cm.triangles = (List*)m->data;
                cn.center = Point();
                // mark1 triangles of component n
                FOREACHVTTRIANGLE(cn.triangles, t, nn) {
                    pointList.appendTail(t->v1());
                    pointList.appendTail(t->v2());
                    pointList.appendTail(t->v3());
                    MARK_BIT(t,1);
                    cn.center += t->getCenter();
                }
                cn.center /= cn.triangles->numels();
                cn.hierarchyTriangles = ot->getTriangleListFromPath(
                        ot->getPathForPointList(pointList, false).toConstPath());
                bool found = false;
                // mark0 triangles in hierarchy of n, that are not part of n
                FOREACHVTTRIANGLE(cn.hierarchyTriangles, t, nn) {
                    if( !IS_BIT(t,1) ) { // t not part of component n
                        MARK_BIT(t,0);
                        found=true;
                    }
                }
                int found_counter = 0;
                if(found) { // other triangles in the hierarchy of n found, overlap possible
                    JMesh::info(" Other triangles in the hierarchy of %d found, overlap possible.\n",i);
                    found = false;
                    cn.inSphereRadius2 = DBL_MAX;
                    cn.outSphereRadius2 = 0;
                    FOREACHVTTRIANGLE(cn.triangles, t, nn) {
                        double d = t->getCenter().squaredDistance(&cn.center);
                        if(d < cn.inSphereRadius2) cn.inSphereRadius2 = d;
                        if(d > cn.outSphereRadius2) cn.outSphereRadius2 = d;
                    }
                    // mark2 all triangles of component n, that contain a point inside m
                    FOREACHVTTRIANGLE(cm.triangles, t, nn) {
                        if( IS_BIT(t,0) && // t was in hierarchy of component n (and could possible overlap with it)
                            ( isPointInComponent(t->v1(), 1, &cn.center, &cn.inSphereRadius2, &cn.outSphereRadius2) ||
                              isPointInComponent(t->v2(), 1, &cn.center, &cn.inSphereRadius2, &cn.outSphereRadius2) ||
                              isPointInComponent(t->v3(), 1, &cn.center, &cn.inSphereRadius2, &cn.outSphereRadius2) )) {
                            MARK_BIT(t,2);
                            need_unlink = true;
                            found = true;
                            found_counter++;
                        }
                    }
                    printf("12\n");
                    if(found) { // triangles of m inside of n found, reverse overlap possible
                        JMesh::info(" %d triangles of %d inside of %d found, searching for triangles of %d in %d.\n", found_counter, j, i, i, j);
                        pointList.removeNodes();
                        // mark3 triangles of m
                        cm.center = Point();
                        FOREACHVTTRIANGLE(cm.triangles, t, nn) {
                            pointList.appendTail(t->v1());
                            pointList.appendTail(t->v2());
                            pointList.appendTail(t->v3());
                            MARK_BIT(t,3);
                            cm.center += t->getCenter();
                        }
                        cm.center /= cm.triangles->numels();
                        cm.inSphereRadius2 = DBL_MAX;
                        cm.outSphereRadius2 = 0;
                        FOREACHVTTRIANGLE(cm.triangles, t, nn) {
                            double d = t->getCenter().squaredDistance(&cm.center);
                            if(d < cm.inSphereRadius2) cm.inSphereRadius2 = d;
                            if(d > cm.outSphereRadius2) cm.outSphereRadius2 = d;
                        }
                        cm.hierarchyTriangles = ot->getTriangleListFromPath(
                                ot->getPathForPointList(pointList, false).toConstPath());
                        found_counter = 0;
                        // mark2 all triangles of component n, that contain a point inside m
                        FOREACHVTTRIANGLE(cn.hierarchyTriangles, t, nn) {
                            if ( IS_BIT(t,1) && // t was part of component n
                                 !IS_BIT(t,2) && // don't recheck if already marked for deletion
                                 ( isPointInComponent((Point*)t->v1(), 3, &cm.center, &cm.inSphereRadius2, &cm.outSphereRadius2) ||
                                   isPointInComponent((Point*)t->v2(), 3, &cm.center, &cm.inSphereRadius2, &cm.outSphereRadius2) ||
                                   isPointInComponent((Point*)t->v3(), 3, &cm.center, &cm.inSphereRadius2, &cm.outSphereRadius2) )) {
                                MARK_BIT(t,2);
                                need_unlink = true;
                                found_counter++;
                            }
                        }
                        JMesh::info(" %d triangles of %d in %d found.\n", found_counter, i, j);
                        cm.hierarchyTriangles->removeNodes();
                        // unmark3 triangles of m again
                        FOREACHVTTRIANGLE(cm.triangles, t, nn) UNMARK_BIT(t,3);
                    } else JMesh::info(" No triangles of %d in %d found.\n", i, j);
                }
                // unmark0 triangles in hierarchy of n again, that are not part of n
                FOREACHVTTRIANGLE(cn.hierarchyTriangles, t, nn) UNMARK_BIT(t,0);
                // mark triangles of n within the minAllowedDistance of any triangle of m for deletion
                double minAllowedDistance2 = minAllowedDistance*minAllowedDistance;
                FOREACHVTTRIANGLE(cm.triangles, t, nn) {
                    if(!IS_BIT(t,2)) { // if the t is not already marked for deletion
                        Point cc = t->getCircleCenter();
                        List *closeTriangles = ot->getTriangleListFromPath(
                                ot->getPathForSphere(cc, minAllowedDistance, false).toConstPath());
                        FOREACHVTTRIANGLE(closeTriangles, t2, mm) {
                            if ( !IS_BIT(t2,2) && // if the t2 is not already marked for deletion
                                 IS_BIT(t2,1) && // is part of component n
                                 ( t2->getCircleCenter().squaredDistance(&cc) <  minAllowedDistance2 )) {
                                // mark both triangles for deletion
                                MARK_BIT(t2,2);
                                need_unlink = true;
                            }
                        }
                        closeTriangles->removeNodes();
                    }
                }
                cn.hierarchyTriangles->removeNodes();
                // delete marked2 triangles, now both components should have at least one boundary loop
                if (need_unlink) {
                    d_boundaries = d_handles = d_shells = 1;
                    int deletion_counter = 0;
                    // remove from component n
                    List *tmp = new List();
                    while( t = (Triangle*) cn.triangles->popHead() ) {
                        if (IS_BIT(t,2)) {
                            this->unlinkTriangle(t);
                            deletion_counter++;
                        } else tmp->appendTail(t);
                    }
                    cn.triangles->joinTailList(tmp);
                    // remove from component m
                    while( t = (Triangle*) cm.triangles->popHead() ) {
                        if (IS_BIT(t,2)) {
                            this->unlinkTriangle(t);
                            deletion_counter++;
                        } else tmp->appendTail(t);
                    }
                    cm.triangles->joinTailList(tmp);
                    JMesh::info(" Deleted %d triangles\n", deletion_counter);
                    ot->removeUnlinkedTriangles();
                    d_boundaries = d_handles = d_shells = 1;
                    removeUnlinkedElements();
                    // check if both components have a boundary loop and join them
                    // TODO: factor 1.5 needed?
//                    JMesh::info(" Joining close boundaries of %d,%d\n", i,j);
//                    if( joinComponentsCloseBoundaries(nl, ml, minAllowedDistance*1.5 ) ) {
//                        // break out of the components loop, reinitialize the components list,
//                        // since the number of the components has changed now
//                        break;
//                    }
                    this->unmarkEverything();
                    this->eulerUpdate();
                    const char *msg = this->checkConnectivity();
                    if(msg)
                        JMesh::error(msg);
                    return 1;
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

bool ExtTriMesh::isPointInComponent(Point *p, char bit, Point *center, double *innerSphereRadius2, double *outerSphereRadius2) {
    double dCenter = p->squaredDistance(center);
    if( dCenter < *innerSphereRadius2 ) return true;
    if( dCenter > *outerSphereRadius2 ) return false;
    TriangleOctree::cursor::path n = ot->getPathForSphere(*p, 1.0, false);
    if ( ot->findClosestPathContainingTriangleMask(n, 1<<bit) ) {
        List *l = ot->getTriangleListFromPath(n.toConstPath());
        double mdist = DBL_MAX;
        bool found = false;
        Triangle *t, *closestTriangle;
        Node *nn;
        FOREACHVTTRIANGLE(l, t, nn) {
            if( (t->mask & 1<<bit) > 0 ) { // part of component
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
