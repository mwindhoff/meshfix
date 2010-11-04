#include "exttrimesh.h"
#include "component.h"
#include "detectIntersections.h"

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

    int ret = this->joinComponentsCloseBoundaries(c1.triangles, c2.triangles, 10);
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
        bool joined = false;
        while (loop2 = (List*) cm.boundaries->popHead()) {
            if(loopsHaveAllVerticesCloserThanDistance(loop1, loop2, joinDistance)
                && closestPair(loop1, loop2, &bv, &bw)
                && joinBoundaryLoops(bv, bw, true, false, false)) {
                ret++;
                joined = true;
                delete(loop2);
                break;
            } else tmp2->appendHead(loop2); // retry this boundary later
        }
        if(!joined) tmp1->appendHead(loop1);
        cm.boundaries->joinTailList(tmp2); // leave boundaries that had no partner
        delete(tmp2);
    }
    // now join all other boundaries
    while ( loop1 = (List*) tmp1->popHead() ) {
        Vertex *bv, *bw;
        while (loop2 = (List*) cm.boundaries->popHead()) {
            if( closestPair(loop1, loop2, &bv, &bw)
                && joinBoundaryLoops(bv, bw, false, false, false)) {
                ret++;
                delete(loop2);
                break;
            }
            delete(loop2);
        }
        delete(loop1);
    }
    this->fillSmallBoundaries(this->E.numels());
    cn.clear();
    cm.clear();
    delete(tmp1);
    return ret;
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

bool ExtTriMesh::decoupleFirstFromSecondComponent(double minAllowedDistance, unsigned max_iterations) {
    bool quiet = JMesh::quiet;
    int iteration_counter = 0;
    short innerBit = 4, outerBit = 5;
    Triangle *t;
    Node *n;
    Vertex *v, *v1, *v2, *v3;
    ExtTriMesh *outer, *inner; // temporary triangulation for intermediate cleaning etc.
    std::map<Vertex*, Point> shift;
    std::map<Vertex*, int> numNormals;
    if(this->shells() != 2) JMesh::error("Must have exactly 2 components.\n");
    outer = (ExtTriMesh*) this->extractFirstShell(); // extract outer (first) component
    inner = (ExtTriMesh*) this->extractFirstShell(); // extract inner (second) component
    this->joinTailTriangulation(outer);
    delete(outer);
    // dilate the inner component by d, to contrain the minAllowedDistance
    inner->dilate(minAllowedDistance);
    JMesh::quiet = true; inner->clean(); JMesh::quiet = quiet;
    while(iteration_counter++ < max_iterations) {
        this->selectAllTriangles(outerBit); // mark outer component
        inner->selectAllTriangles(innerBit); // mark inner component
        this->joinHeadTriangulation(inner);
        JMesh::info("Iteration %d\n", iteration_counter);
        // mark0 triangles of the outer component which are inside of the inner component
        if(!this->markTrianglesInsideComponent(0, outerBit, innerBit)) break; // finished
        // we have triangles of the outer component that are inside or too close to the inner component
        inner = (ExtTriMesh*) this->extractFirstShell(); // extract inner component
        shift.clear();
        numNormals.clear();
        FOREACHVERTEX(v, n) { // initialize shift with (0,0,0)
            shift.insert(std::pair<Vertex*, Point>(v,Point()));
            numNormals.insert(std::pair<Vertex*, int>(v,0));
        }
        FOREACHTRIANGLE(t, n) if(IS_BIT(t,0)) { // compute shift for affected vertices
            UNMARK_BIT(t,0);
            Point normal = t->getNormal();
            v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
            shift[v1] += normal; shift[v2] += normal; shift[v3] += normal;
            numNormals[v1]++; numNormals[v2]++; numNormals[v3]++;
        }
        // compute shift as mean normal of surrounding triangles
        for(std::map<Vertex*, Point>::iterator it = shift.begin(); it != shift.end(); ++it) {
            v = it->first;
            if(int num = numNormals[v]) *v += it->second/num;
        }
        this->unmarkEverything();
        JMesh::report_progress("Cleaning ...");
        JMesh::quiet = true;
        this->clean(); // and clean and repair it (because the shift could have produced new intersections ...)
        this->checkAndRepair();
        JMesh::quiet = quiet;
        JMesh::report_progress("");
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
