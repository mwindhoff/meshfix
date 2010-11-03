#include "exttrimesh.h"
#include "component.h"
#include "detectIntersections.h"

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

/* Assumes that the Triangulation consists of exactly 2 components, each having no selfintersections.
   If they overlap, they will be joined and the overlapping parts will be deleted. */
int ExtTriMesh::joinOverlappingComponentPair() {
    this->deselectTriangles();
    List *components = this->getComponentsList();
    if( components->numels() > 2 ) JMesh::error("Triangulation consists of more than 2 components.\n");
    if( components->numels() < 2 ) { JMesh::info("Only 1 component, nothing joined.\n"); return 0; }
    // select the intersecting triangles, which form the boundaries
    if(!this->selectIntersectingTriangles()) {
        JMesh::info("Component pair doesn't overlap. Nothing joined.\n");
        return 0;
    }
    // Identify the two most distant triangles (which are assumed to be not part of the overlap),
    // to determine which chunks to keep after removing the intersecting triangles
    ComponentStruct c1((List*) components->popHead()), c2((List*) components->popHead());
    Triangle *t1, *t2;
    Node *n;
    FOREACHVTTRIANGLE(c1.triangles, t1, n) if(IS_VISITED(t1)) break;
    FOREACHVTTRIANGLE(c2.triangles, t2, n) if(IS_VISITED(t2)) break;
    c1.vertices = c1.getVertices(2);
    c2.vertices = c2.getVertices(2);
    Vertex *v1, *v2;
    this->mostDistantPartner(t1->v1(), c2.vertices, &v2);
    this->mostDistantPartner(t2->v1(), c1.vertices, &v1);
    c1.clear();
    c2.clear();
    t1 = v1->e0->t1;
    t2 = v2->e0->t1;
    // now delete the boundaries and have at least 2 shells with boundaries
    this->removeSelectedTriangles();
    if(!t1 || !t2) JMesh::error("Algorithm using most distant points didn't work ...\n");
    this->selectConnectedComponent(t1, false);
    this->selectConnectedComponent(t2, false);
    this->invertSelection();
    this->removeSelectedTriangles();
    // remove more triangles close to the overlap, to make the holes bigger for better joining
    this->selectBoundaryTriangles();
    this->growSelection();
    this->removeSelectedTriangles();
    this->removeSmallestComponents(2);
    c1 = ComponentStruct(t1);
    c2 = ComponentStruct(t2);
    int ret = this->joinComponentsCloseBoundaries(c1.triangles, c2.triangles, DBL_MAX);
    this->eulerUpdate();
    this->fillSmallBoundaries(0, true, false);
    return ret;
}

/* Assumes that the Triangulation consists of exactly 2 components, each having no selfintersections.
   If they overlap, they will be joined and the overlapping parts will be deleted. */
int ExtTriMesh::joinOverlappingComponentPair2() {
    this->deselectTriangles();
    List *components = this->getComponentsList();
    if( components->numels() > 2 ) JMesh::error("Triangulation consists of more than 2 components.\n");
    if( components->numels() < 2 ) { JMesh::info("Only 1 component, nothing joined.\n"); return 0; }
    // select the intersecting triangles, which form the boundaries
    List *first = (List*) components->head()->data;
    Triangle *t; Node *n;
    // mark4 first component
    FOREACHVTTRIANGLE(first, t, n) {
        MARK_BIT(t,4);
        MARK_BIT(t->v1(),4);
        MARK_BIT(t->v2(),4);
        MARK_BIT(t->v3(),4);
    }
    delete(first);
    delete(components);
    // mark5 second component
    FOREACHTRIANGLE(t, n) if(!IS_BIT(t,4)) {
        MARK_BIT(t,5);
        MARK_BIT(t->v1(),5);
        MARK_BIT(t->v2(),5);
        MARK_BIT(t->v3(),5);
    }
    this->markTrianglesInsideComponent(6, 5, 4);
    this->markTrianglesInsideComponent(7, 4, 5);
    FOREACHTRIANGLE(t, n) if(t->mask & (1<<6 | 1<<7)) t->mask = 1;
    this->removeSelectedTriangles();
    this->eulerUpdate();
    if(this->shells() != 2) return 0;
    FOREACHTRIANGLE(t, n) if(!IS_BIT(t,5)) break;
    ComponentStruct c1(t);
    FOREACHTRIANGLE(t, n) if(!IS_BIT(t,4)) break;
    ComponentStruct c2(t);

    int ret = this->joinComponentsCloseBoundaries(c1.triangles, c2.triangles, 3);
    if(!ret) {
        this->selectBoundaryTriangles();
        this->removeSelectedTriangles();
        FOREACHTRIANGLE(t, n) if(!IS_BIT(t,5)) break;
        c1 = ComponentStruct(t);
        FOREACHTRIANGLE(t, n) if(!IS_BIT(t,4)) break;
        c2 = ComponentStruct(t);
        ret = this->joinComponentsCloseBoundaries(c1.triangles, c2.triangles, 10);
        if(!ret) return 0;
    }
    this->eulerUpdate();
    this->fillSmallBoundaries(0, true, false);
    this->unmarkEverything();
    return ret;
}

int ExtTriMesh::joinComponentsCloseBoundaries(List *nl, List *ml, double joinDistance) {
    ComponentStruct cn(nl), cm(ml);
    cn.initializeBoundaries();
    cm.initializeBoundaries();
    List *loop1, *loop2;
    int ret = 0;
    List *tmp1 = new List();
    // first join all boundaries that are closer than the joinDistance
    while ( loop1 = (List*) cn.boundaries->popHead() ) {
        List *tmp2 = new List();
        Vertex *bv, *bw;
        while (loop2 = (List*) cm.boundaries->popHead()) {
            if(loopsHaveAllVerticesCloserThanDistance(loop1, loop2, joinDistance)
                && closestPair(loop1, loop2, &bv, &bw)
                && joinBoundaryLoops(bv, bw, false, true, false)) {
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

double ExtTriMesh::mostDistantPartner(Vertex *v, List *l, Vertex **distantPartner) {
    Node *m;
    Vertex *w;
    double adist, maxdist = 0;
    FOREACHVVVERTEX(l, w, m) {
        if ((adist = w->squaredDistance(v)) > maxdist) {
            maxdist=adist;
            *distantPartner = w;
        }
    }
    return maxdist;
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

int ExtTriMesh::moveVerticesInwards(Point &componentCenter, std::map<Vertex *, Point> &origin, double stepsize, double distance) {
    List todo(new di_cell(this)), cells;
    di_cell *c, *c2;
    int i = 0;
    int ret = 0;
    double stepsize2 = stepsize*stepsize;
    while (c = (di_cell *)todo.popHead()) {
        if (i > DI_MAX_NUMBER_OF_CELLS || c->triangles.numels() <= 10 || (c->Mp-c->mp).length() < distance) cells.appendHead(c);
        else {
            i++;
            JMesh::report_progress(NULL);
            c2 = c->fork();
            if (!c->containsBothShells(1,2)) delete(c); else todo.appendTail(c);
            if (!c2->containsBothShells(1,2)) delete(c2); else todo.appendTail(c2);
        }
    }
    JMesh::report_progress("");
    Node *n, *m;
    Triangle *t, *t2;
    Vertex *v1, *v2;
    double distance2 = distance*distance;
    while (c = (di_cell *)cells.popHead()) {
        std::set<Vertex*> vertices;
        FOREACHVTTRIANGLE((&c->triangles), t, n) { vertices.insert(t->v1()); vertices.insert(t->v2()); vertices.insert(t->v3()); }
        for(std::set<Vertex*>::const_iterator itv1 = vertices.begin(); itv1 != vertices.end(); ++itv1) {
            v1 = *itv1;
            if(IS_BIT(v1, 2)) {
                FOREACHVTTRIANGLE((&c->triangles), t, n) if(IS_BIT(t, 1)) {
                    Point center = t->getCircleCenter();
                    double radius2 = MAX(center.squaredDistance(t->v1()),distance2);
                    if(center.squaredDistance(v1) < radius2) {
                        Point n = t->getNormal()*distance, p = center + n;
                        v2->intersectionWithPlane(v1, &origin[v1], &p, &n);
                        if(v1->squaredDistance(&origin[v1]) > v2->squaredDistance(&origin[v1]))
                            *v1 += *v2 - *v1;
                        UNMARK_BIT(v1,2);
                    }
                }
            }
        }
    }
    FOREACHVERTEX(v1, n) if(IS_BIT(v1,2)) {
        double dist = origin[v1].distance(v1);
        if(dist > stepsize) *v1 += (origin[v1]-(*v1))*(stepsize/dist);
        else { *v1 = origin[v1]; UNMARK_BIT(v1,2); }
        ret++;
    }
    return ret;
}

bool ExtTriMesh::decoupleSecondFromFirstComponent(double minAllowedDistance, unsigned max_iterations) {
    int iteration_counter = 0;
    // save origin of each point
    Vertex *v, *v1, *v2, *v3;
    ExtTriMesh *tmptin; // temporary triangulation for intermediate cleaning
    std::map<Vertex*, Point> origin, shift;
    std::map<Vertex*, int> numNormals;
    Node *n;
    tmptin = (ExtTriMesh*) this->extractShell((Triangle*) T.head()->data);
    // dilate the inner component by d, to maintain the minAllowedDistance
    tmptin->dilate(minAllowedDistance);
    this->joinTailTriangulation(tmptin);
    T.appendTail(T.popHead());
    delete(tmptin);
    FOREACHVERTEX(v, n)
        origin.insert(std::pair<Vertex*, Point>(v,Point(*v)));
    while(iteration_counter++ < max_iterations) {
        // switch the order of the components (by moving one triangle to the head of the list)
        T.appendHead(T.popTail());
        List *components = this->getComponentsList();
        if(components->numels() != 2) JMesh::error("Must have exactly 2 components.\n");
        List *first = (List*) components->head()->data;
        Triangle *t;
        // mark4 first (inner) component
        FOREACHVTTRIANGLE(first, t, n) {
            MARK_BIT(t,4);
            MARK_BIT(t->v1(),4);
            MARK_BIT(t->v2(),4);
            MARK_BIT(t->v3(),4);
        }
        delete(first);
        delete(components);
        // mark5 second (outer) component
        FOREACHTRIANGLE(t, n) if(!IS_BIT(t,4)) {
            MARK_BIT(t,5);
            MARK_BIT(t->v1(),5);
            MARK_BIT(t->v2(),5);
            MARK_BIT(t->v3(),5);
        }
        shift.clear();
        numNormals.clear();
        FOREACHVERTEX(v, n) { // initialize shift with (0,0,0)
            shift.insert(std::pair<Vertex*, Point>(v,Point()));
            numNormals.insert(std::pair<Vertex*, int>(v,0));
        }
        JMesh::info("Iteration %d\n", iteration_counter);
        // mark0 triangles of componenent2 (mark5) which are inside of component1 (mark4)
        if(!this->markTrianglesInsideComponent(0, 5, 4)) {
            // mark0 vertices of componenent2 (mark5) which are inside of component1 (mark4)
//            if(!this->moveTooCloseVerticesOutwards(minAllowedDistance, 5, 4)) {
                // delete first shell (since we didn't change it we are only interested in the second)
                FOREACHTRIANGLE(t, n) if(IS_BIT(t,4)) { // get first shell
                    this->unmarkEverything();
                    this->removeShell(t);
                    this->eulerUpdate();
                    break;
                }
                break; // finished
//            }
        } else { // we have triangles of the component2 that are inside or too close to the component1
            FOREACHTRIANGLE(t, n) {
                if(IS_BIT(t,0) && IS_BIT(t,5)) { // compute shift for affected vertices
                    Point normal = t->getNormal();
                    v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
                    shift[v1] += normal; shift[v2] += normal; shift[v3] += normal;
                    numNormals[v1]++; numNormals[v2]++; numNormals[v3]++;
                }
                UNMARK_BIT(t,0);
            }
            // compute shift as mean normal of surrounding triangles
            FOREACHVERTEX(v, n) if(int num = numNormals[v]) *v += shift[v]/num;
        }
        FOREACHTRIANGLE(t, n) if(IS_BIT(t,5)) { // get first triangle of component1
            this->unmarkEverything();
            tmptin = (ExtTriMesh*) this->extractShell(t); // extract this component
            JMesh::quiet = true;
            tmptin->clean(); // and clean and repair it (because the shift could have produced new intersections ...)
            tmptin->checkAndRepair();
            JMesh::quiet = false;
            this->deselectTriangles();
            this->joinTailTriangulation(tmptin); // insert it again
            delete(tmptin);
            this->eulerUpdate();
            JMesh::quiet = true; this->checkAndRepair(); JMesh::quiet = false;
            break;
        }
    }
    if(iteration_counter < max_iterations) return true;
    return false;
}

int ExtTriMesh::markTrianglesInsideComponent(short insideMarkBit, short componentMarkBit1, short componentMarkBit2) {
    di_cell *c = new di_cell(this), *c2, *tmp;
    Triangle *t;
    Node *n;
    List todo(c), cells, tmptl;
    // keep only triangles of the two components
    while(t = (Triangle*) c->triangles.popHead())
        if(t->mask & (1<<componentMarkBit1 | 1<<componentMarkBit2))
            tmptl.appendHead(t);
    c->triangles.joinTailList(&tmptl);
    int ncells = 0;
    int ret = 0;
    // get smallest cells containing at least both shells
    while (c = (di_cell *)todo.popHead()) {
        if (ncells > 10*DI_MAX_NUMBER_OF_CELLS || c->triangles.numels() <= 100) cells.appendHead(c);
        else {
            JMesh::report_progress(NULL);
            tmp = new di_cell(*c);
            c2 = c->fork();
            if (!c->containsBothShells(componentMarkBit1, componentMarkBit2) ||
                !c2->containsBothShells(componentMarkBit1, componentMarkBit2)) {
                delete(c);
                delete(c2);
                cells.appendHead(tmp);
            } else {
                ncells++;
                todo.appendTail(c);
                todo.appendTail(c2);
                delete(tmp);
            }
        }
    }
    JMesh::report_progress("");
    // if no intersections, then all triangles
    int nintersections = this->selectIntersectingTriangles();
    std::set<Vertex*> vertices1, vertices2;
    short outsideMarkBit = 0;
    // get first unused bit
    unsigned char mask = 1<<componentMarkBit1 | 1<<componentMarkBit2 | 1<<insideMarkBit | 1;
    while(1<<outsideMarkBit & mask) outsideMarkBit++;
    unsigned char decidedMask = 1<<insideMarkBit | 1<<outsideMarkBit;
    while(c = (di_cell*) cells.popHead()) {
        JMesh::report_progress("%d%%", (int)round(((double)(ncells - cells.numels()))/ncells*100));
        // get vertices of triangles of component1 in the cell
        Vertex *vt;
        FOREACHVTTRIANGLE((&c->triangles), t, n) {
            bool comp1 = IS_BIT(t, componentMarkBit1);
            vt = t->v1(); if(!comp1) vertices2.insert(vt); else if(!(vt->mask & decidedMask)) vertices1.insert(vt);
            vt = t->v2(); if(!comp1) vertices2.insert(vt); else if(!(vt->mask & decidedMask)) vertices1.insert(vt);
            vt = t->v3(); if(!comp1) vertices2.insert(vt); else if(!(vt->mask & decidedMask)) vertices1.insert(vt);
        }
        // decide for each vertex whether inside or outside
        for(std::set<Vertex*>::const_iterator i = vertices1.begin(); i != vertices1.end(); ++i) {
            Vertex *v = *i, *w, *closest;
            // search for yet undecided (unmarked) connected regions
            if(!(v->mask & decidedMask)) {
                // get closest vertex of the other component
                double d, dmin = DBL_MAX;
                for(std::set<Vertex*>::const_iterator j = vertices2.begin(); j != vertices2.end(); ++j) {
                    w = *j;
                    if(!IS_VISITED(w)) {
                        MARK_VISIT(w);
                        d = v->squaredDistance(w);
                        if(d < dmin) { dmin = d; closest = w; }
                    }
                }
                for(std::set<Vertex*>::const_iterator j = vertices2.begin(); j != vertices2.end(); ++j)
                    UNMARK_VISIT(*j);
                if(dmin == DBL_MAX) { /*v->printPoint(); */continue; }
                // get the mean normal at the closest vertex
                Point meanNormal;
                FOREACHVTTRIANGLE(closest->VT(), t, n) meanNormal += t->getNormal();
                // decide whether it is inside or outside (using the mean normal plane)
                bool isInside = meanNormal.squaredLength() ? meanNormal*(*v - *closest) < 0 : false;
                // find an unselected triangle of having vertex v
                List *l = v->VT();
                while(t = (Triangle*) l->popHead()) if(!IS_VISITED(t)) { todo.appendHead(t); break; };
                delete(l);
                if(!todo.numels()) continue; // no unselected triangle in neighborhood
                // spread the decision to all connected triangles, stop at selected triangles (== intersections)
                short markbit = isInside ? insideMarkBit : outsideMarkBit;
                while(t = (Triangle*) todo.popHead()) {
                    MARK_BIT(t, markbit);
                    MARK_BIT(t->v1(), markbit); MARK_BIT(t->v2(), markbit); MARK_BIT(t->v3(), markbit);
                    Triangle *t1 = t->t1(), *t2 = t->t2(), *t3 = t->t3();
                    // stop at selected triangles (== markbit 0 == intersections)
                    if(t1 && !(t1->mask & (decidedMask | 1<<0))) todo.appendHead(t1);
                    if(t2 && !(t2->mask & (decidedMask | 1<<0))) todo.appendHead(t2);
                    if(t3 && !(t3->mask & (decidedMask | 1<<0))) todo.appendHead(t3);
                }
            }
        }
        delete(c);
        vertices1.clear();
    }
    JMesh::report_progress("");
    FOREACHTRIANGLE(t, n) {
        if(IS_BIT(t, componentMarkBit1)) {
            if(IS_VISITED(t) || IS_BIT(t, insideMarkBit)) {
                UNMARK_VISIT(t);
                MARK_BIT(t, insideMarkBit);
                ret++;
            } else UNMARK_BIT(t, outsideMarkBit);
            t->v1()->mask = 0; t->v2()->mask = 0; t->v3()->mask = 0;
        } else UNMARK_VISIT(t);
    }
    JMesh::info("Number of triangles inside: %d\n", ret);
    return ret;
}

int ExtTriMesh::moveTooCloseVerticesOutwards(double minAllowedDistance, short componentMarkBit1, short componentMarkBit2) {
    di_cell *c = new di_cell(this), *c2;
    Triangle *t, *t2;
    Node *n, *m;
    List todo(c), cells, tmptl;
    // keep only triangles of the two components
    while(t = (Triangle*) c->triangles.popHead())
        if(t->mask & (1<<componentMarkBit1 | 1<<componentMarkBit2))
            tmptl.appendHead(t);
    c->triangles.joinTailList(&tmptl);
    int ncells = 0;
    int ret = 0;
    // cellsize = sqrt(d^2+d^2+d^2) = sqrt(3*d^3)
    double cellsize2 = 4*3*minAllowedDistance*minAllowedDistance;
    // get smallest cells containing at least both shells
    while (c = (di_cell *)todo.popHead()) {
        JMesh::report_progress(NULL);
        if (ncells > DI_MAX_NUMBER_OF_CELLS || c->triangles.numels() <= 10 || (c->Mp-c->mp).squaredLength() < cellsize2 )
            cells.appendHead(c);
        else {
            ncells++;
            JMesh::report_progress(NULL);
            c2 = c->fork();
            if (c->containsBothShells(componentMarkBit1, componentMarkBit2))
                todo.appendTail(c);
            else delete(c);
            if (c2->containsBothShells(componentMarkBit1, componentMarkBit2))
                todo.appendTail(c2);
            else delete(c2);
        }
    }
    double minAllowedDistance2 = minAllowedDistance*minAllowedDistance;
    std::set<Vertex *> vertices;
    std::map<Vertex *, Point> shift;
    std::map<Vertex *, double> minDist2;
    Vertex *v;
    FOREACHVERTEX(v, n) {
        minDist2[v] = minAllowedDistance2;
        shift[v] = Point();
    }
    while(c = (di_cell*) cells.popHead()) {
        vertices.clear();
        JMesh::report_progress(NULL);
        FOREACHVTTRIANGLE((&c->triangles), t, n) if(IS_BIT(t, componentMarkBit1)) {
            vertices.insert(t->v1());
            vertices.insert(t->v2());
            vertices.insert(t->v3());
        }
        while(t = (Triangle*) c->triangles.popHead()) if(IS_BIT(t, componentMarkBit2)) {
            for(std::set<Vertex *>::iterator i = vertices.begin(); i != vertices.end(); ++i) {
                double dist2 = t->pointTriangleSquaredDistance(*i);
                if (dist2 < minDist2[*i]) {
                    minDist2[*i] = dist2;
                    MARK_VISIT(*i);
                }
            }
        }
        delete(c);
    }
    FOREACHVERTEX(v, n) if(IS_VISITED(v)) {
        List *vtl = v->VT();
        FOREACHVTTRIANGLE(vtl, t2, m) shift[v] += t2->getNormal();
        shift[v] /= vtl->numels();
        delete(vtl);
    }
    FOREACHVERTEX(v, n) if(IS_VISITED(v)) {
        ret++;
        *v += shift[v]*MAX((minAllowedDistance-sqrt(minDist2[v])), MAX(0.1*minAllowedDistance, 1));
        UNMARK_VISIT(v);
    }
    JMesh::report_progress("");
    JMesh::info("Number of too close vertices: %d\n", ret);
    return ret;
}

void ExtTriMesh::dilate(double d) {
    Vertex *v; Node *n, *m; Triangle *t;
    std::map<Vertex *, Point> shift;
    int nsteps = MAX((int) d,1);
    double step = d/(double)nsteps;
    for(int i = 0; i < nsteps; i++) {
        FOREACHVERTEX(v, n)  {
            shift[v] = Point();
            List *vtl = v->VT();
            FOREACHVTTRIANGLE(vtl, t, m) shift[v] += t->getNormal();
            shift[v] /= vtl->numels()/step;
            delete(vtl);
        }
        FOREACHVERTEX(v, n) *v += shift[v];
        this->clean();
    }
}
