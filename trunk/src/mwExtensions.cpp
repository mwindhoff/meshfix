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
}

struct ComponentStruct {
    List *triangles;
    List *hierarchyTriangles;
    Point center;
    double inSphereRadius2;
    double outSphereRadius2;
    void initializeRadiiAndCenter() {
        Triangle *t; Node *n;
        inSphereRadius2 = DBL_MAX;
        outSphereRadius2 = 0;
        center = Point();
        FOREACHVTTRIANGLE(triangles, t, n) center += t->getCenter();
        center /= triangles->numels();
        FOREACHVTTRIANGLE(triangles, t, n) {
            double d1 = t->v1()->squaredDistance(&center);
            double d2 = t->v2()->squaredDistance(&center);
            double d3 = t->v3()->squaredDistance(&center);
            inSphereRadius2 = MIN(inSphereRadius2, MIN(d1, MIN(d2, d3)));
            outSphereRadius2 = MAX(outSphereRadius2, MAX(d1, MAX(d2, d3)));
        }
    }
    void markBit(unsigned b) {
        Triangle *t; Node *n;
        FOREACHVTTRIANGLE(triangles, t, n) MARK_BIT(t,b);
    }
    void unmarkBit(unsigned b) {
        Triangle *t; Node *n;
        FOREACHVTTRIANGLE(triangles, t, n) UNMARK_BIT(t,b);
    }
};
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
    ComponentStruct cn, cm; // component structs
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
                cn.initializeRadiiAndCenter();
                // mark1 triangles of component n
                FOREACHVTTRIANGLE(cn.triangles, t, nn) {
                    pointList.appendTail(t->v1());
                    pointList.appendTail(t->v2());
                    pointList.appendTail(t->v3());
                    MARK_BIT(t,1);
                }
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
                    if(found) { // triangles of m inside of n found, reverse overlap possible
                        JMesh::info(" %d triangles of %d inside of %d found, searching for triangles of %d in %d.\n", found_counter, j, i, i, j);
                        pointList.removeNodes();
                        // mark3 triangles of m
                        FOREACHVTTRIANGLE(cm.triangles, t, nn) {
                            pointList.appendTail(t->v1());
                            pointList.appendTail(t->v2());
                            pointList.appendTail(t->v3());
                            MARK_BIT(t,3);
                        }
                        cm.initializeRadiiAndCenter();
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
                        cm.unmarkBit(3);
                    } else JMesh::info(" No triangles of %d in %d found.\n", i, j);
                }
                pointList.removeNodes();
                // unmark0 triangles in hierarchy of n again, that are not part of n
                FOREACHVTTRIANGLE(cn.hierarchyTriangles, t, nn) UNMARK_BIT(t,0);
                // mark for deletion triangles within the minAllowedDistance of any not yet deleted triangle
                double minAllowedDistance2 = minAllowedDistance*minAllowedDistance;
                FOREACHVTTRIANGLE(cm.triangles, t, mm) {
                    if(!IS_BIT(t,2)) { // if the t is not already marked for deletion
                        Point cc = t->getCircleCenter();
                        FOREACHVTTRIANGLE(cn.triangles, t2, nn) {
                            if ( !IS_BIT(t2,2) && // if the t2 is not already marked for deletion
                                 (!IS_BIT(t,3) || !IS_BIT(t2,3)) && // one of them is not yet marked for later deletion
                                 ( t2->getCircleCenter().squaredDistance(&cc) <  minAllowedDistance2 )) {
                                // mark3 both triangles (for later deletion)
                                MARK_BIT(t,3);
                                MARK_BIT(t2,3);
                                need_unlink = true;
                            }
                        }
                    }
                }
                // mark2 all triangles for deletion which are marked3
                FOREACHVTTRIANGLE(cn.triangles, t, nn) if(IS_BIT(t,3)) MARK_BIT(t,2);
                FOREACHVTTRIANGLE(cm.triangles, t, nn) if(IS_BIT(t,3)) MARK_BIT(t,2);
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
                    this->removeUnlinkedElements();
                    this->unmarkEverything();
                    this->eulerUpdate();
                    return 1;
                    if( joinComponentsCloseBoundaries(cn.triangles, cm.triangles, minAllowedDistance) )
                        printf("true\n");
                    const char *msg = this->checkConnectivity();
                    if(msg)
                        JMesh::error(msg);
                    // check if both components have a boundary loop and join them
                    // TODO: factor 1.5 needed?
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
    Node *n, *m, *nn;
    Vertex *v,*w, *v1, *v2, *v3;
    Triangle *t, *t1, *t2, *t3;
    List *vertexList = new List(), *nBoundaryLoops = new List(), *mBoundaryLoops = new List(), *one_loop = new List();
    List **nBoundaryLoopsArray;
    // create list of vertices of component n
    t = (Triangle*)nl->head()->data;
    List todo(t); MARK_BIT(t,2);
    while (todo.numels()) {
        t = (Triangle *)todo.popHead();
        t1 = t->t1(); t2 = t->t2(); t3 = t->t3();
        v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
        if (!IS_VISITED(v1)) {MARK_VISIT(v1); vertexList->appendHead(v1);}
        if (!IS_VISITED(v2)) {MARK_VISIT(v2); vertexList->appendHead(v2);}
        if (!IS_VISITED(v3)) {MARK_VISIT(v3); vertexList->appendHead(v3);}

        if (t1 != NULL && !IS_BIT(t1,2)) {MARK_BIT(t1,2); todo.appendHead(t1);}
        if (t2 != NULL && !IS_BIT(t2,2)) {MARK_BIT(t2,2); todo.appendHead(t2);}
        if (t3 != NULL && !IS_BIT(t3,2)) {MARK_BIT(t3,2); todo.appendHead(t3);}
    }
    FOREACHVTTRIANGLE(nl, t, n) UNMARK_BIT(t,2);
    printf("vln: %d, n: %d\n", vertexList->numels(), nl->numels());
    // create list of list boundary loops of component n (= list of vertices)
    FOREACHVVVERTEX(vertexList, v, n) {
        // find next vertex of an unmarked boundary
        if (!IS_BIT(v,2) && v->isOnBoundary()) {
            w = v;
            do { // mark all vertices at this boundary
                one_loop->appendHead(w);
                MARK_BIT(w,2);
                w = w->nextOnBoundary();
            } while (w != v);
            nBoundaryLoops->joinTailList(one_loop);
        }
    }
    FOREACHVVVERTEX(vertexList, v, n) {UNMARK_BIT(v,1); UNMARK_BIT(v,2);}
    vertexList->removeNodes();
    // create list of vertices of component m
    t = (Triangle*)ml->head()->data;
    todo.appendHead(t); MARK_BIT(t,2);
    while (todo.numels()) {
        t = (Triangle *)todo.popHead();
        t1 = t->t1(); t2 = t->t2(); t3 = t->t3();
        v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
        if (!IS_VISITED(v1)) {MARK_VISIT(v1); vertexList->appendHead(v1);}
        if (!IS_VISITED(v2)) {MARK_VISIT(v2); vertexList->appendHead(v2);}
        if (!IS_VISITED(v3)) {MARK_VISIT(v3); vertexList->appendHead(v3);}

        if (t1 != NULL && !IS_BIT(t1,2)) {MARK_BIT(t1,2); todo.appendHead(t1);}
        if (t2 != NULL && !IS_BIT(t2,2)) {MARK_BIT(t2,2); todo.appendHead(t2);}
        if (t3 != NULL && !IS_BIT(t3,2)) {MARK_BIT(t3,2); todo.appendHead(t3);}
    }
    FOREACHVTTRIANGLE(ml, t, n) UNMARK_BIT(t,2);
    printf("vln: %d, n: %d\n", vertexList->numels(), ml->numels());
    // create list of list boundary loops of component m (= list of vertices)
    FOREACHVVVERTEX(vertexList, v, n) {
        // find next vertex of an unmarked boundary
        if (!IS_BIT(v,2) && v->isOnBoundary()) {
            w = v;
            do { // mark all vertices at this boundary
                one_loop->appendHead(w);
                MARK_BIT(w,2);
                w = w->nextOnBoundary();
            } while (w != v);
            mBoundaryLoops->joinTailList(one_loop);
        }
    }
    FOREACHVVVERTEX(vertexList, v, n) {UNMARK_BIT(v,1); UNMARK_BIT(v,2);}
    vertexList->removeNodes();
    nBoundaryLoopsArray = (List **)nBoundaryLoops->toArray();
    unsigned nNumLoops = nBoundaryLoops->numels();
    bool allVerticesHaveClosePartner = false;
    double maxDistanceToJoin2 = maxDistanceToJoin*maxDistanceToJoin;
    bool ret = false;
    int j = 0;
    printf("nNumLoops %d, mNumLoops %d\n", nNumLoops, mBoundaryLoops->numels());
    for (unsigned i = 0; i < nNumLoops; i++) {
        printf("nLoops %d\n", i);
        nn = mBoundaryLoops->head();
        do {
            printf("mLoops %d\n", j++);
            allVerticesHaveClosePartner = true;
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
            JMesh::info(" allVerticesHaveClosePartner: %d\n", allVerticesHaveClosePartner);
            if ( allVerticesHaveClosePartner ) {
                joinBoundaryLoops(v, w, true, true, true);
                ret = true;
                mBoundaryLoops->freeCell(nn);
                break;
            }
            nn = nn->next();
        } while ( nn != NULL );
    }
    delete(nBoundaryLoopsArray);
    while ((one_loop=(List *)nBoundaryLoops->popHead())!=NULL) delete one_loop;
    while ((one_loop=(List *)mBoundaryLoops->popHead())!=NULL) delete one_loop;
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

