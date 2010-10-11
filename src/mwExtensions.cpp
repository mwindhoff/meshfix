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
    ot->writeTrianglesToLeafs();
}


/**
 * Joins close or overlapping components of triangles.
 * Algorithm sketch:
 * assume all triangles unmarked
 * for all pairs of components (m,n):
 *  mark1 all triangles of n
 *  mark0 all triangles that are not in n and that are in the hierarchy of n
 *  mark2 all triangles of m that are marked0 and that have a point inside of n
 *  mark2 all triangles that are in the hierarchy of m, have a point inside of m and are marked1 (==part of n)
 *  for all triangles on boundary of m
 *   mark2 all triangles of n within the minAllowedDistance
 *  delete marked2 triangles
 *  get the list of boundary loops of m and n
 *   find a boundary within the minimal joining distance
 *    find closes pair of vertices and join them with filling, refining and fairing.
 */
int ExtTriMesh::joinCloseOrOverlappingComponents( double minAllowedDistance ) {
    JMesh::info("Joining components that overlap or are closer than %g:\n", minAllowedDistance);
    Node *n, *m, *nn, *mm, *nnn;
    Triangle *t, *t2, *t3;
    ExtTriMesh *two;
    unsigned i = 0, j = 0;
    // fill components list
    bool need_unlink = false, joined = false;
    do {
        ComponentStruct *cn, *cm; // component structs
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
            cn = new ComponentStruct((Triangle*)two->T.head()->data);
            cn->markBit(1); // mark1 triangles of component n
            t = (Triangle*)((List*)m->data)->head()->data;
            two->append(new ExtTriMesh(t));
            this->removeShell(t);
            cm = new ComponentStruct();
            cm->triangles = new List();
            FOREACHVTTRIANGLE((&two->T), t, nn) if (!IS_BIT(t,1)) cm->triangles->appendHead(t);
            two->initializeOctree();
            cn->initializeRadiiAndCenter();
            cn->vertices = cn->getVertices(2);
            cn->hierarchyTriangles = new List();
            two->ot->fillTriangleListFromPath(cn->hierarchyTriangles,
                    two->ot->getPathForPointList(*cn->vertices, false).toConstPath());
            int found_counter = 0;
            // mark0 triangles in hierarchy of n, that are not part of n
            FOREACHVTTRIANGLE(cn->hierarchyTriangles, t, nn) {
                if( !IS_BIT(t,1) ) { // t not part of component n
                    MARK_BIT(t,0);
                    found_counter++;
                }
            }
            if(found_counter) { // other triangles in the hierarchy of n found, overlap possible
                JMesh::info(" (%d,%d): Some triangles could possibly overlap.\n",i+1,j+1, found_counter);
                found_counter = 0;
                // mark2 all triangles of component m, that contain a point inside n
                FOREACHVTTRIANGLE(cm->triangles, t, nn) {
                    if( IS_BIT(t,0) && // t was in hierarchy of component n (and could possible overlap with it)
                        ( two->isPointInComponent(t->v1(), cn) ||
                          two->isPointInComponent(t->v2(), cn) ||
                          two->isPointInComponent(t->v3(), cn) )) {
                        MARK_BIT(t,2);
                        need_unlink = true;
                        found_counter++;
                    }
                }
                if(found_counter) { // triangles of m inside of n found, reverse overlap possible
                    JMesh::info(" (%d,%d): %d overlapping triangles of component %d.\n", i+1, j+1, found_counter, j+1);
                    // mark3 triangles of m
                    cm->vertices = cm->getVertices(3);
                    cm->initializeRadiiAndCenter();
                    cm->hierarchyTriangles = new List();
                    TriangleOctree::cursor::path p = two->ot->getPathForPointList(*cm->vertices, false);
                    two->ot->fillTriangleListFromPath(cm->hierarchyTriangles, p.toConstPath());
                    found_counter = 0;
                    // mark2 all triangles of component n, that contain a point inside m
                    FOREACHVTTRIANGLE(cm->hierarchyTriangles, t, nn) {
                        if ( IS_BIT(t,1) && // t was part of component n
                             !IS_BIT(t,2) && // don't recheck if already marked for deletion
                             ( two->isPointInComponent(t->v1(), cm) ||
                               two->isPointInComponent(t->v2(), cm) ||
                               two->isPointInComponent(t->v3(), cm) )) {
                            MARK_BIT(t,2);
                            need_unlink = true;
                            found_counter++;
                        }
                    }
                    JMesh::info(" (%d,%d): %d overlapping triangles of component %d.\n", i+1, j+1, found_counter, i+1);
                    cm->hierarchyTriangles->removeNodes();
                } else JMesh::info(" (%d,%d): No overlap.\n", i+1, j+1);
            }
            // unmark0 triangles in hierarchy of n again
            FOREACHVTTRIANGLE(cn->hierarchyTriangles, t, nn) UNMARK_BIT(t,0);
            // mark for deletion triangles within the minAllowedDistance of any not yet deleted triangle
            double minAllowedDistance2 = minAllowedDistance*minAllowedDistance;
            FOREACHVTTRIANGLE(cm->triangles, t, mm) {
                if(!IS_BIT(t,2)) { // if the t is not already marked for deletion
                    Point cc = t->getCircleCenter();
                    FOREACHVTTRIANGLE(cn->triangles, t2, nn) {
                        // mark3 all connected triangles (for later deletion)
                        if ( t2->v1()->squaredDistance(&cc) <  minAllowedDistance2 ) {
                            List *vt = t2->v1()->VT();
                            FOREACHVTTRIANGLE(vt, t3, nnn) MARK_BIT(t3,3);
                            delete(vt);
                            MARK_BIT(t,3);
                            need_unlink = true;
                        } else if ( t2->v2()->squaredDistance(&cc) <  minAllowedDistance2 ) {
                            List *vt = t2->v2()->VT();
                            FOREACHVTTRIANGLE(vt, t3, nnn) MARK_BIT(t3,3);
                            delete(vt);
                            MARK_BIT(t,3);
                            need_unlink = true;
                        } else if ( t2->v3()->squaredDistance(&cc) <  minAllowedDistance2 ) {
                            List *vt = t2->v3()->VT();
                            FOREACHVTTRIANGLE(vt, t3, nnn) MARK_BIT(t3,3);
                            delete(vt);
                            MARK_BIT(t,3);
                            need_unlink = true;
                        }
                        // mark3 all connected triangles (for later deletion)
                        if (IS_BIT(t,3)) {
                            List *vt = t->v1()->VT();
                            FOREACHVTTRIANGLE(vt, t3, nnn) MARK_BIT(t3,4);
                            vt->removeNodes();
                            vt = t->v2()->VT();
                            FOREACHVTTRIANGLE(vt, t3, nnn) MARK_BIT(t3,4);
                            vt->removeNodes();
                            vt = t->v3()->VT();
                            FOREACHVTTRIANGLE(vt, t3, nnn) MARK_BIT(t3,4);
                            delete(vt);
                        }
                    }
                }
            }
            found_counter = 0;
            // mark2 all triangles for deletion which are marked3
            FOREACHVTTRIANGLE(cn->triangles, t, nn) {
                if(!IS_BIT(t,2) && IS_BIT(t,3)) {
                    MARK_BIT(t,2);
                    found_counter++;
                }
            }
            FOREACHVTTRIANGLE(cm->triangles, t, nn) {
                if(!IS_BIT(t,2) && (IS_BIT(t,3) || IS_BIT(t,4))) {
                    MARK_BIT(t,2);
                    found_counter++;
                }
            }
            if( found_counter )
                JMesh::info(" (%d,%d) %d triangles are closer than %g.\n", i+1, j+1, found_counter, minAllowedDistance);
            // delete marked2 triangles, now both components should have at least one boundary loop
            if (need_unlink) {
                d_boundaries = d_handles = d_shells = 1;
                int deletion_counter = 0;
                // remove from component n
                List *tmp = new List();
                while( t = (Triangle*) cn->triangles->popHead() ) {
                    if (IS_BIT(t,2)) {
                        two->unlinkTriangle(t);
                        deletion_counter++;
                    } else tmp->appendTail(t);
                }
                cn->triangles->joinTailList(tmp);
                JMesh::info(" Deleted %d triangles of component %d\n", deletion_counter, i+1);
                deletion_counter = 0;
                // remove from component m
                while( t = (Triangle*) cm->triangles->popHead() ) {
                    if (IS_BIT(t,2)) {
                        two->unlinkTriangle(t);
                        deletion_counter++;
                    } else tmp->appendTail(t);
                }
                cm->triangles->joinTailList(tmp);
                JMesh::info(" Deleted %d triangles of component %d\n", deletion_counter, j+1);
                two->ot->removeUnlinkedTriangles();
                two->d_boundaries = two->d_handles = two->d_shells = 1;
                two->removeUnlinkedElements();
                two->unmarkEverything();
                // check if both components have a boundary loop and join them
                if( two->joinComponentsCloseBoundaries(cn->triangles, cm->triangles, 2*minAllowedDistance) ) {
                    this->append(two);
                    delete(two);
                    this->eulerUpdate();
                    joined = true;
                    JMesh::info(" (%d,%d): Components successfully joined.\n", i+1, j+1);
                    i=0;
                    j=0;
                    cm->clear();
                    cn->clear();
                    break;
                }
            }
            this->append(two);
            delete(two);
            this->eulerUpdate();
            if( !joined ) {
                j++;
                break;
            }
        }
        cm->clear();
        cn->clear();
        // delete components list
        FOREACHNODE(*components, n) delete((List *)n->data);
    } while (j < this->n_shells); // until there were no more triangles to delete
    this->unmarkEverything();
    return 0;
}

int ExtTriMesh::joinComponentsCloseBoundaries(List *nl, List *ml, double joinDistance) {
    ComponentStruct cn(nl), cm(ml);
    cn.initializeBoundaries();
    cm.initializeBoundaries();
    List *loop1, *loop2;
    bool ret = false;
    List *tmp1 = new List();
    // first join all boundaries that are closer than the joinDistance
    while ( loop1 = (List*) cn.boundaries->popHead() ) {
        List *tmp2 = new List();
        Vertex *bv, *bw;
        while (loop2 = (List*) cm.boundaries->popHead()) {
            if(loopsHaveAllVerticesCloserThanDistance(loop1, loop2, joinDistance)
                && closestPair(loop1, loop2, &bv, &bw)
                && joinBoundaryLoops(bv, bw, false, true, true)) {
                ret++;
                loop2->removeNodes();
                break;
            } else tmp2->appendHead(loop2); // retry this boundary later
        }
        tmp1->appendHead(loop1);
        cm.boundaries->joinTailList(tmp2); // leave boundaries that had no partner
    }
    // now join all other boundaries
    while ( loop1 = (List*) tmp1->popHead() ) {
        Vertex *bv, *bw;
        while (loop2 = (List*) cm.boundaries->popHead()) {
            if( closestPair(loop1, loop2, &bv, &bw)
                && joinBoundaryLoops(bv, bw, false, true, true)) {
                ret++;
                loop2->removeNodes();
                break;
            }
        }
        loop1->removeNodes();
    }
    cn.clear();
    cm.clear();
    return ret;
}

bool ExtTriMesh::isPointInComponent(Vertex *v, ComponentStruct *c) {
    double dCenter = v->squaredDistance(&c->center);
    if( dCenter < c->inSphereRadius2 ) return true;
    if( dCenter > c->outSphereRadius2 ) return false;
    Vertex *w;
    this->getClosestPartner(v, c->vertices, &w);
    Point n = w->e0->t1->getNormal();
    if (n.isNull()) return false;
    double d = n*(*v) - n*(*w);
    return d < 0;
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
