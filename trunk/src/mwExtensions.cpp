#include "exttrimesh.h"
#include "component.h"
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
//    ot->writeTrianglesToLeafs();
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
    JMesh::info("Joining components that overlap or are closer than %g:\n", minAllowedDistance);
    Node *n, *m, *nn, *mm, *nnn;
    Triangle *t, *t2, *t3;
    ExtTriMesh *two;
    unsigned i = 0, j = 0;
    // fill components list
    bool need_unlink = false, joined = false;
    ComponentStruct cn, cm; // component structs
    do {
        need_unlink = false;
        List *components = getComponentsList();
        JMesh::info(" Number of components: %d\n", components->numels() );
        // handle all possible pairs of components
        unsigned ii = 0, jj = 0;
        FOREACHNODE(*components, n) {
            if( ii < i ) continue; // skip already checked components
            m=n->next();
            while( jj++ < j && m ) m=m->next();
            if (!m) { // no component left
                i++;
                j=i+1;
                break;
            }
            joined = false;
            JMesh::info(" Checking component pair (%d,%d) for overlap.\n", i+1, j+1);
            // create triangulation of both components
            t = (Triangle*)((List*)n->data)->head()->data;
            two = new ExtTriMesh(t);
            this->removeShell(t);
            cn = ComponentStruct((Triangle*)two->T.head()->data);
            cn.markBit(1); // mark1 triangles of component n
            List pointList(two->V); // list of vertices of cn
            t = (Triangle*)((List*)m->data)->head()->data;
            two->append(new ExtTriMesh(t));
            this->removeShell(t);
            cm.triangles = new List();
            FOREACHVTTRIANGLE((&two->T), t, nn) if (!IS_BIT(t,1)) cm.triangles->appendHead(t);
            two->initializeOctree();
            cn.initializeRadiiAndCenter();
            cn.hierarchyTriangles = two->ot->getTriangleListFromPath(
                    two->ot->getPathForPointList(pointList, false).toConstPath());
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
                JMesh::info(" (%d,%d): Possible overlap.\n",i+1,j+1);
                found = false;
                // mark2 all triangles of component n, that contain a point inside m
                FOREACHVTTRIANGLE(cm.triangles, t, nn) {
                    if( IS_BIT(t,0) && // t was in hierarchy of component n (and could possible overlap with it)
                        ( two->isPointInComponent(t->v1(), 1, &cn.center, &cn.inSphereRadius2, &cn.outSphereRadius2) ||
                          two->isPointInComponent(t->v2(), 1, &cn.center, &cn.inSphereRadius2, &cn.outSphereRadius2) ||
                          two->isPointInComponent(t->v3(), 1, &cn.center, &cn.inSphereRadius2, &cn.outSphereRadius2) )) {
                        MARK_BIT(t,2);
                        need_unlink = true;
                        found = true;
                        found_counter++;
                    }
                }
                if(found) { // triangles of m inside of n found, reverse overlap possible
                    JMesh::info(" (%d,%d): %d overlapping triangles of component %d.\n", i+1, j+1, found_counter, j+1);
                    pointList.removeNodes();
                    // mark3 triangles of m
                    pointList = *cm.getVertices();
                    cm.initializeRadiiAndCenter();
                    cm.hierarchyTriangles = two->ot->getTriangleListFromPath(
                            two->ot->getPathForPointList(pointList, false).toConstPath());
                    found_counter = 0;
                    // mark2 all triangles of component n, that contain a point inside m
                    FOREACHVTTRIANGLE(cn.hierarchyTriangles, t, nn) {
                        if ( IS_BIT(t,1) && // t was part of component n
                             !IS_BIT(t,2) && // don't recheck if already marked for deletion
                             ( two->isPointInComponent((Point*)t->v1(), 3,
                                                       &cm.center, &cm.inSphereRadius2, &cm.outSphereRadius2) ||
                               two->isPointInComponent((Point*)t->v2(), 3,
                                                       &cm.center, &cm.inSphereRadius2, &cm.outSphereRadius2) ||
                               two->isPointInComponent((Point*)t->v3(), 3,
                                                       &cm.center, &cm.inSphereRadius2, &cm.outSphereRadius2) )) {
                            MARK_BIT(t,2);
                            need_unlink = true;
                            found_counter++;
                        }
                    }
                    JMesh::info(" (%d,%d): %d overlapping triangles of component %d.\n", i+1, j+1, found_counter, i+1);
                    cm.hierarchyTriangles->removeNodes();
                    // unmark3 triangles of m again
                    cm.unmarkBit(3);
                } else JMesh::info(" (%d,%d): No overlap.\n", i+1, j+1);
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
                        if ( true ) { // one of them is not yet marked for later deletion
                            // mark3 all connected triangles (for later deletion)
                            if ( t2->v1()->squaredDistance(&cc) <  minAllowedDistance2 ) {
                                FOREACHVTTRIANGLE(t2->v1()->VT(), t3, nnn) MARK_BIT(t3,3);
                                MARK_BIT(t,3);
                                need_unlink = true;
                            }
                            if ( t2->v2()->squaredDistance(&cc) <  minAllowedDistance2 ) {
                                FOREACHVTTRIANGLE(t2->v2()->VT(), t3, nnn) MARK_BIT(t3,3);
                                MARK_BIT(t,3);
                                need_unlink = true;
                            }
                            if ( t2->v3()->squaredDistance(&cc) <  minAllowedDistance2 ) {
                                FOREACHVTTRIANGLE(t2->v3()->VT(), t3, nnn) MARK_BIT(t3,3);
                                MARK_BIT(t,3);
                                need_unlink = true;
                            }
                            // mark3 all connected triangles (for later deletion)
                            if (IS_BIT(t,3)) {
                                FOREACHVTTRIANGLE(t->v1()->VT(), t3, nnn) MARK_BIT(t3,4);
                                FOREACHVTTRIANGLE(t->v2()->VT(), t3, nnn) MARK_BIT(t3,4);
                                FOREACHVTTRIANGLE(t->v3()->VT(), t3, nnn) MARK_BIT(t3,4);
                            }
                        }
                    }
                }
            }
            // mark2 all triangles for deletion which are marked3
            FOREACHVTTRIANGLE(cn.triangles, t, nn) if(IS_BIT(t,3)) MARK_BIT(t,2);
            FOREACHVTTRIANGLE(cm.triangles, t, nn) if(IS_BIT(t,3) || IS_BIT(t,4)) MARK_BIT(t,2);
            cn.hierarchyTriangles->removeNodes();
            // delete marked2 triangles, now both components should have at least one boundary loop
            if (need_unlink) {
                d_boundaries = d_handles = d_shells = 1;
                int deletion_counter = 0;
                // remove from component n
                List *tmp = new List();
                while( t = (Triangle*) cn.triangles->popHead() ) {
                    if (IS_BIT(t,2)) {
                        two->unlinkTriangle(t);
                        deletion_counter++;
                    } else tmp->appendTail(t);
                }
                cn.triangles->joinTailList(tmp);
                JMesh::info(" Deleted %d triangles of component %d\n", deletion_counter, i+1);
                deletion_counter = 0;
                // remove from component m
                while( t = (Triangle*) cm.triangles->popHead() ) {
                    if (IS_BIT(t,2)) {
                        two->unlinkTriangle(t);
                        deletion_counter++;
                    } else tmp->appendTail(t);
                }
                cm.triangles->joinTailList(tmp);
                JMesh::info(" Deleted %d triangles of component %d\n", deletion_counter, j+1);
                two->ot->removeUnlinkedTriangles();
                two->d_boundaries = two->d_handles = two->d_shells = 1;
                two->removeUnlinkedElements();
                two->unmarkEverything();
                // check if both components have a boundary loop and join them
                if( two->joinComponentsCloseBoundaries(cn.triangles, cm.triangles, 2*minAllowedDistance) ) {
                    this->append(two);
                    this->eulerUpdate();
                    joined = true;
                    JMesh::info(" (%d,%d): Components successfully joined.\n", i+1, j+1);
                    i=0;
                    j=0;
                    break;
                }
            }
            this->append(two);
            this->eulerUpdate();
            if( !joined ) {
                j++;
                break;
            }
        }
        // delete components list
        FOREACHNODE(*components, n) delete((List *)n->data);
    } while (j < this->n_shells); // until there were no more triangles to delete
    this->unmarkEverything();
    return 0;
}

int ExtTriMesh::joinComponentsCloseBoundaries(List *nl, List *ml, double maxDistanceToJoin) {
    ComponentStruct cn(nl), cm(ml);
    cn.initializeBoundaries();
    cm.initializeBoundaries();
    List *loop, *loop2;
    bool ret = false;
    while ( loop = (List*) cn.boundaries->popHead() ) {
        List *tmp = new List();
        Vertex *bv, *bw;
        while (loop2 = (List*) cm.boundaries->popHead()) {
            if(loopsHaveAllVerticesCloserThanDistance(loop, loop2, maxDistanceToJoin)) {
                if( closestPair(loop, loop2, &bv, &bw) ) {
                    if(joinBoundaryLoops(bv, bw, false, true, true)) {
                        ret++;
                        break;
                    }
                }
                loop2->removeNodes();
                break;
            } else tmp->appendHead(loop2); // retry this boundary later
        }
        loop->removeNodes();
        cm.boundaries->joinTailList(tmp); // leave boundaries that had no partner
    }
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

bool ExtTriMesh::loopsHaveAllVerticesCloserThanDistance(List *loop, List *loop2, const double &distance) {
    Node *n, *m;
    Vertex *v, *w;
    const double d2 = distance*distance;
    if (loop->numels() && loop2->numels()) {
        FOREACHVVVERTEX(loop, v, n) {
            bool foundClosePartner = false;
            FOREACHVVVERTEX(loop2, w, m) {
                if( v->squaredDistance(w) < d2 ) {
                    foundClosePartner = true;
                    break;
                }
            }
            if (!foundClosePartner) return false;
        }
        return true;
    }
    return false;
}

double ExtTriMesh::closestPair(List *l1, List *l2, Vertex **closest1, Vertex **closest2)
{
    Node *n, *m;
    Vertex *v,*w;
    double adist, mindist = DBL_MAX;
    FOREACHVVVERTEX(l1, v, n) {
        FOREACHVVVERTEX(l2, w, m) {
            if ((adist = w->squaredDistance(v)) < mindist) {
                mindist=adist;
                *closest1 = v;
                *closest2 = w;
            }
        }
    }
    return mindist;
}

double ExtTriMesh::getClosestPartner(Vertex *v, List *l, Vertex **closestParnter) {
    Node *m;
    Vertex *w;
    double adist, mindist = DBL_MAX;
    FOREACHVVVERTEX(l, w, m) {
        if ((adist = w->squaredDistance(v)) < mindist) {
            mindist=adist;
            *closestParnter = w;
        }
    }
    return mindist;
}
