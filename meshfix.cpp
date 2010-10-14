#include "exttrimesh.h"
#include <string.h>
#include <stdlib.h>
#include "jrs_predicates.h"

const char *input_filename;

// Simulates the ASCII rounding error
void ExtTriMesh::asciiAlign()
{
 char outname[2048];
 Vertex *v;
 Node *m;
 float a;
 FOREACHVERTEX(v, m)
 {
  sprintf(outname,"%f",v->x); sscanf(outname,"%f",&a); v->x = a;
  sprintf(outname,"%f",v->y); sscanf(outname,"%f",&a); v->y = a;
  sprintf(outname,"%f",v->z); sscanf(outname,"%f",&a); v->z = a;
 }
}


// Return TRUE if the triangle is exactly degenerate
inline bool isDegenerateEdge(Edge *e)
{
 return ((*(e->v1))==(*(e->v2)));
}

bool isDegenerateTriangle(Triangle *t)
{
 double xy1[2], xy2[2], xy3[2];
 xy1[0] = t->v1()->x; xy1[1] = t->v1()->y;
 xy2[0] = t->v2()->x; xy2[1] = t->v2()->y;
 xy3[0] = t->v3()->x; xy3[1] = t->v3()->y;
 if (orient2d(xy1, xy2, xy3)!=0.0) return false;
 xy1[0] = t->v1()->y; xy1[1] = t->v1()->z;
 xy2[0] = t->v2()->y; xy2[1] = t->v2()->z;
 xy3[0] = t->v3()->y; xy3[1] = t->v3()->z;
 if (orient2d(xy1, xy2, xy3)!=0.0) return false;
 xy1[0] = t->v1()->z; xy1[1] = t->v1()->x;
 xy2[0] = t->v2()->z; xy2[1] = t->v2()->x;
 xy3[0] = t->v3()->z; xy3[1] = t->v3()->x;
 if (orient2d(xy1, xy2, xy3)!=0.0) return false;
 return true;
}

Edge *getLongestEdge(Triangle *t)
{
 double l1 = t->e1->squaredLength();
 double l2 = t->e2->squaredLength();
 double l3 = t->e3->squaredLength();
 if (l1>=l2 && l1>=l3) return t->e1;
 if (l2>=l1 && l2>=l3) return t->e2;
 return t->e3;
}

// Iterate on all the selected triangles as long as possible.
// Keep the selection only on the degeneracies that could not be removed.
// Return the number of degeneracies that could not be removed
int ExtTriMesh::swapAndCollapse()
{
 Node *n;
 Triangle *t;

 if (epsilon_angle != 0.0)
 {
  FOREACHTRIANGLE(t, n) UNMARK_VISIT(t);
  JMesh::quiet = true; removeDegenerateTriangles(); JMesh::quiet = false;
  int failed = 0;
  FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) failed++;
  return failed;
 }

 List triangles;
 Edge *e;
 const int MAX_ATTEMPTS = 10;

 FOREACHTRIANGLE(t, n) t->info=0;

 // VISIT2 means that the triangle is in the list
 FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
 {
  UNMARK_VISIT(t);
  if (isDegenerateTriangle(t)) {triangles.appendTail(t); MARK_VISIT2(t);}
 }

 while ((t=(Triangle *)triangles.popHead())!=NULL)
 {
  UNMARK_VISIT2(t);
  if (t->isLinked())
  {
   if (isDegenerateEdge(t->e1)) t->e1->collapse();
   else if (isDegenerateEdge(t->e2)) t->e2->collapse();
   else if (isDegenerateEdge(t->e3)) t->e3->collapse();
   else if ((e=getLongestEdge(t))!=NULL)
   {
    if (e->swap())
    {
     t=e->t1;
     if (isDegenerateTriangle(t) && !IS_VISITED2(t) && ((long int)t->info < MAX_ATTEMPTS))
     {triangles.appendTail(t); MARK_VISIT2(t); t->info = (void *)(((long int)t->info)+1);}
     t=e->t2;
     if (isDegenerateTriangle(t) && !IS_VISITED2(t) && ((long int)t->info < MAX_ATTEMPTS))
     {triangles.appendTail(t); MARK_VISIT2(t); t->info = (void *)(((long int)t->info)+1);}
    }
   }
  }
 }

 removeUnlinkedElements();

 int failed=0;
 // This should check only on actually processed triangles
 FOREACHTRIANGLE(t, n) if (isDegenerateTriangle(t)) {failed++; MARK_VISIT(t);}

 JMesh::info("%d degeneracies selected\n",failed);
 return failed;
}

// returns true on success

bool ExtTriMesh::removeDegenerateTriangles(int max_iters, int num_to_keep)
{
 int n, iter_count = 0;

 printf("Removing degeneracies...\n");
 while ((++iter_count) <= max_iters && swapAndCollapse())
 {
  for (n=1; n<iter_count; n++) growSelection();
  removeSelectedTriangles();
  removeSmallestComponents(num_to_keep);
  JMesh::quiet = true; fillSmallBoundaries(E.numels()); JMesh::quiet = false;
  asciiAlign();
 }

 if (iter_count > max_iters) return false;
 return true;
}

bool appendCubeToList(Triangle *t0, List& l)
{
 if (!IS_VISITED(t0) || IS_VISITED2(t0)) return false;

 Triangle *t, *s;
 Vertex *v;
 List triList(t0);
 MARK_VISIT2(t0);
 double minx=DBL_MAX, maxx=-DBL_MAX, miny=DBL_MAX, maxy=-DBL_MAX, minz=DBL_MAX, maxz=-DBL_MAX;

 while(triList.numels())
 {
  t = (Triangle *)triList.popHead();
  v = t->v1();
  minx=MIN(minx,v->x); miny=MIN(miny,v->y); minz=MIN(minz,v->z);
  maxx=MAX(maxx,v->x); maxy=MAX(maxy,v->y); maxz=MAX(maxz,v->z);
  v = t->v2();
  minx=MIN(minx,v->x); miny=MIN(miny,v->y); minz=MIN(minz,v->z);
  maxx=MAX(maxx,v->x); maxy=MAX(maxy,v->y); maxz=MAX(maxz,v->z);
  v = t->v3();
  minx=MIN(minx,v->x); miny=MIN(miny,v->y); minz=MIN(minz,v->z);
  maxx=MAX(maxx,v->x); maxy=MAX(maxy,v->y); maxz=MAX(maxz,v->z);
  if ((s = t->t1()) != NULL && !IS_VISITED2(s) && IS_VISITED(s)) {triList.appendHead(s); MARK_VISIT2(s);}
  if ((s = t->t2()) != NULL && !IS_VISITED2(s) && IS_VISITED(s)) {triList.appendHead(s); MARK_VISIT2(s);}
  if ((s = t->t3()) != NULL && !IS_VISITED2(s) && IS_VISITED(s)) {triList.appendHead(s); MARK_VISIT2(s);}
 }

 l.appendTail(new Point(minx, miny, minz));
 l.appendTail(new Point(maxx, maxy, maxz));
 return true;
}

bool isVertexInCube(Vertex *v, List& loc)
{
 Node *n;
 Point *p1, *p2;
 FOREACHNODE(loc, n)
 {
  p1 = (Point *)n->data; n=n->next(); p2 = (Point *)n->data;
  if (!(v->x < p1->x || v->y < p1->y || v->z < p1->z ||
      v->x > p2->x || v->y > p2->y || v->z > p2->z)) return true;
 }

 return false;
}

void ExtTriMesh::selectTrianglesInCubes()
{
 Triangle *t;
 Vertex *v;
 Node *n;
 List loc;
 FOREACHTRIANGLE(t, n) appendCubeToList(t, loc);
 FOREACHVERTEX(v, n) if (isVertexInCube(v, loc)) MARK_VISIT(v);
 FOREACHTRIANGLE(t, n)
 {
  UNMARK_VISIT2(t);
  if (IS_VISITED(t->v1()) || IS_VISITED(t->v2()) || IS_VISITED(t->v3())) MARK_VISIT(t);
 }
 FOREACHVERTEX(v, n) UNMARK_VISIT(v);
 loc.freeNodes();
}

// returns true on success

bool ExtTriMesh::removeSelfIntersections(int max_iters, int number_components_to_keep)
{
 int n, iter_count = 0;

 printf("Removing self-intersections...\n");
 while ((++iter_count) <= max_iters && selectIntersectingTriangles())
 {
  for (n=1; n<iter_count; n++) growSelection();
  removeSelectedTriangles();
  removeSmallestComponents(number_components_to_keep);
  JMesh::quiet = true; fillSmallBoundaries(E.numels()); JMesh::quiet = false;
  asciiAlign();
  selectTrianglesInCubes();
 }

 if (iter_count > max_iters) return false;
 return true;
}

bool ExtTriMesh::removeSelfIntersections2(int max_iterations, int number_components_to_keep)
{
    int iteration_counter = 0, remove_and_fill_counter = 0, smooth_counter = 0, grow_counter = 0;
    printf("Removing self-intersections (using advanced method)...\n");
    int nintersecting = 0, nintersecting_new = 0;
    deselectTriangles();
    invertSelection();
    JMesh::info("Stage: Remove and Fill (1)\n");
    while (true)
    {
        iteration_counter++;
        asciiAlign();
        if((nintersecting_new = selectIntersectingTriangles(10)) > 0) {
            remove_and_fill_counter++;
            // remove intersecting triangles
            removeSelectedTriangles();
            // remove smallest shells
            removeSmallestComponents(number_components_to_keep);
            // fill, refine, fair, keep new triangles selected
            JMesh::quiet=true; fillSmallBoundaries(E.numels(), true); JMesh::quiet=false;
            // grow selection, recheck selection for intersections
            growSelection();
            if (nintersecting != nintersecting_new && remove_and_fill_counter < max_iterations*2) {
                // the last iteration resulted in different holes as before
                nintersecting = nintersecting_new;
                continue;
            }
        } else {
            deselectTriangles();
            if(!selectIntersectingTriangles())
                return true; // we have reached the end
            continue;
        }
        remove_and_fill_counter = 0;
        JMesh::info("Stage: Laplacian Smooth (%d)\n", smooth_counter+1);
        // next step is smoothing
        deselectTriangles();
        removeSmallestComponents(number_components_to_keep);
        if(!selectIntersectingTriangles()) continue;
        if(smooth_counter++ < max_iterations) {
            JMesh::info("Laplacian smoothing of selected triangles.\n");
            // increase region to smooth
            for( int i = 0; i < smooth_counter; i++) growSelection();
            // smooth with 1 step, keep selection
            JMesh::quiet=true; laplacianSmooth(); JMesh::quiet=false;
            growSelection();
            nintersecting = 0;
            JMesh::info("Stage: Remove and Fill (%d)\n", iteration_counter+1);
            continue;
        }
        smooth_counter = 0;
        JMesh::info("Stage: Grow selection, Remove and Fill (%d)\n", grow_counter+1);
        deselectTriangles();
        removeSmallestComponents(number_components_to_keep);
        if(selectIntersectingTriangles()) {
            for (int i=0; i < grow_counter+1; i++)
                growSelection();
            removeSelectedTriangles();
            removeSmallestComponents(number_components_to_keep);
            JMesh::quiet=true; fillSmallBoundaries(E.numels(), true); JMesh::quiet=false;
            if (++grow_counter >= max_iterations) break;
            JMesh::info("Stage: Remove and Fill (%d)\n", iteration_counter+1);
        }
    }
    return false;
}


bool ExtTriMesh::isDegeneracyFree()
{
 Node *n;
 Triangle *t;

 if (epsilon_angle != 0.0)
 {FOREACHTRIANGLE(t, n) if (t->isDegenerate()) return false;}
 else
 {FOREACHTRIANGLE(t, n) if (isDegenerateTriangle(t)) return false;}

 return true;
}


// returns true on success

bool ExtTriMesh::clean(int max_iters, int inner_loops, int number_components_to_keep)
{
 bool ni, nd;

 deselectTriangles();
 invertSelection();

 for (int n=0; n<max_iters; n++)
 {
  printf("********* ITERATION %d *********\n",n);
  nd=removeDegenerateTriangles(inner_loops);
  deselectTriangles(); invertSelection();
  ni=removeSelfIntersections2(inner_loops, number_components_to_keep);
  if (ni && nd && isDegeneracyFree()) return true;
 }

 return false;
}


double closestPair(List *bl1, List *bl2, Vertex **closest_on_bl1, Vertex **closest_on_bl2)
{
 Node *n, *m;
 Vertex *v,*w;
 double adist, mindist = DBL_MAX;

 FOREACHVVVERTEX(bl1, v, n)
  FOREACHVVVERTEX(bl2, w, m)
   if ((adist = w->squaredDistance(v))<mindist)
   {
    mindist=adist;
    *closest_on_bl1 = v;
    *closest_on_bl2 = w;
   }

 return mindist;
}

/**
 * Joins the closest components, that have boundaries (holes).
 */
bool joinClosestComponents(ExtTriMesh *tin, bool justconnect = false, bool refine = true, bool fair = true) {
    Vertex *v,*w, *gv, *gw;
    Triangle *t, *s;
    Node *n;
    List triList, boundary_loops, *one_loop;
    List **bloops_array;
    int i, j, numloops;

    i=0;
    // delete info of all triangles
    FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;
    // initialize info of all triangles with their component number starting by 1.
    FOREACHVTTRIANGLE((&(tin->T)), t, n) {
        if (t->info == NULL) {
            i++;
            triList.appendHead(t);
            t->info = (void *)i;
            while(triList.numels()) {
                t = (Triangle *)triList.popHead();
                if ((s = t->t1()) != NULL && s->info == NULL) {triList.appendHead(s); s->info = (void *)i;}
                if ((s = t->t2()) != NULL && s->info == NULL) {triList.appendHead(s); s->info = (void *)i;}
                if ((s = t->t3()) != NULL && s->info == NULL) {triList.appendHead(s); s->info = (void *)i;}
            }
        }
    }
    // if less then 2 components
    if (i<2) {
        // unset info again
        FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;
        JMesh::info("Mesh is a single component. Nothing done.");
        return false;
    }
    // copy triangle component number to the vertices
    FOREACHVTTRIANGLE((&(tin->T)), t, n) {
        t->v1()->info = t->v2()->info = t->v3()->info = t->info;
    }
    // create list boundary loop lists (= lists of connected vertices on a boundary)
    FOREACHVVVERTEX((&(tin->V)), v, n) {
        // find next vertex of an unmarked boundary
        if (!IS_VISITED2(v) && v->isOnBoundary()) {
            w = v;
            one_loop = new List;
            // mark all vertices at this boundary
            do {
                one_loop->appendHead(w);
                MARK_VISIT2(w);
                w = w->nextOnBoundary();
            } while (w != v);
            boundary_loops.appendHead(one_loop);
        }
    }
    FOREACHVVVERTEX((&(tin->V)), v, n) UNMARK_VISIT2(v);

    bloops_array = (List **)boundary_loops.toArray();
    numloops = boundary_loops.numels();

    int numtris = tin->T.numels();
    double adist, mindist=DBL_MAX;

    gv=NULL;
    for (i=0; i<numloops; i++) {
        for (j=0; j<numloops; j++) {
            // if i,j are indices of vertices of different boundary loops, search for the closes pair of vertices and update mindist
            if (((Vertex *)bloops_array[i]->head()->data)->info != ((Vertex *)bloops_array[j]->head()->data)->info) {
                adist = closestPair(bloops_array[i], bloops_array[j], &v, &w);
                if (adist<mindist) {mindist=adist; gv=v; gw=w;}
            }
        }
    }
    if (gv!=NULL) tin->joinBoundaryLoops(gv, gw, justconnect, refine, fair);

    FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;
    FOREACHVVVERTEX((&(tin->V)), v, n) v->info = NULL;

    free(bloops_array);
    while ((one_loop=(List *)boundary_loops.popHead())!=NULL) delete one_loop;

    return (gv!=NULL);
}

//#define DISCLAIMER

void usage()
{
 printf("%s v%s - by %s.\n=====================================================\n", JMesh::app_name, JMesh::app_version, JMesh::app_authors);
 printf("USAGE: meshfix <file1> [<file2>] [OPTIONS]\n");
 printf("  Processes <meshfile1> and saves the result to <file1>_fixed.off.\n");
 printf("  An optionally passed <file2> is merged with the first one.\n");
 printf("OPTIONS:\n");
 printf(" -a <epsilon_angle>  Allowed range: 0 < epsilon_angle < 2, default: 0 (degrees).\n");
 printf(" -n <n>              Only the <n> biggest input components are kept.\n");
 printf(" -j                  Join 2 biggest components if they overlap, remove the rest.\n");
// printf(" -j <d>              Join components closer than <d> or overlapping.\n");
 printf(" -jc                 Join the closest pair of components.\n");
 printf(" -u <steps>          Uniform remeshing of the whole mesh, steps > 0\n");
 printf(" --decouple          Treat 1st file as inner, 2nd file as outer component.\n");
 printf("                     Resolve intersections by moving outer triangles outward.\n");
 printf(" --no-clean          Don't clean.\n");
 printf(" -w                  Result is saved in VRML1.0 format instead of OFF.\n");
 printf(" -s                  Result is saved in STL     format instead of OFF.\n");
 printf(" -o <output>         Set the output filename (without extension).\n");
 printf("Accepted input formats are OFF, PLY and STL.\n  Other formats are supported only partially.\n");
 printf("See http://jmeshlib.sourceforge.net for details on supported formats.\n");
 printf("\nIf MeshFix is used for research purposes, please cite the following paper:\n");
 printf("\n   M. Attene.\n   A lightweight approach to repairing digitized polygon meshes.\n   The Visual Computer, 2010. (c) Springer.\n");
 printf("\nHIT ENTER TO EXIT.\n");
 getchar();
 exit(0);
}

char *createFilename(const char *iname, const char *subext, const char *newextension)
{
 static char tname[2048];
 char *oname = (char *)malloc(strlen(iname)+strlen(subext)+strlen(newextension)+1);
 strcpy(tname, iname);
 for (int n=strlen(tname)-1; n>0; n--) if (tname[n]=='.') {tname[n] = '\0'; break;}
 sprintf(oname,"%s%s%s",tname,subext,newextension);
 return oname;
}

int main(int argc, char *argv[])
{
 char subext[128]="_fixed";
 JMesh::init();
 JMesh::app_name = "MeshFix";
 JMesh::app_version = "1.1-alpha";
 JMesh::app_year = "2010";
 JMesh::app_authors = "Marco Attene, Mirko Windhoff";
 JMesh::app_maillist = "attene@ge.imati.cnr.it, mirko.windhoff@tuebingen.mpg.de";

 ExtTriMesh tin;

#ifdef DISCLAIMER
 printf("\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n");
 printf("This software can be used ONLY with an explicit authorization of the author.\n");
 printf("If you do not have such an authorization, you must delete this software.\n");
 printf("In no event this version of MeshFix can be redistributed.\n");
 printf("\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n");
#endif

 if (argc < 2) usage();

 float par = 0;
 unsigned numberComponentsToKeep = 1;
 bool joinOverlappingComponents = false;
// float minAllowedDistance = 0;
 bool haveJoinClosestComponents = false;
 int uniformRemeshSteps = 0;
 bool haveDecouple = false;
 bool clean = true;
 bool save_vrml = false;
 bool save_stl = false;
 bool haveOutputFile = false;
 const char *outputFile;
 for (int i=2; i<argc; i++)
 {
  if (!strcmp(argv[i], "-a"))
  {
   if (i<argc-1) par = (float)atof(argv[i+1]); else par = 0;
   if (par < 0) JMesh::error("Epsilon angle must be > 0.\n");
   if (par > 2) JMesh::error("Epsilon angle must be < 2 degrees.\n");
   tin.epsilon_angle = par;
   if (tin.epsilon_angle)
   {
    JMesh::acos_tolerance = asin((M_PI*tin.epsilon_angle)/180.0);
	printf("Fixing asin tolerance to %e\n",JMesh::acos_tolerance);
    i++;
   }
  }
  else if (!strcmp(argv[i], "-n")) {
      if (i<argc-1) {
          numberComponentsToKeep = atoi(argv[i+1]);
          if (numberComponentsToKeep < 1)
              JMesh::error("# components to keep must be >= 1.\n");
          else
              i++;
      }
  }
  else if (!strcmp(argv[i], "-w")) save_vrml = true;
  else if (!strcmp(argv[i], "-s")) save_stl = true;
  else if (!strcmp(argv[i], "-j")) joinOverlappingComponents = true; /*{
      if (i<argc-1) {
          minAllowedDistance = atof(argv[i+1]);
          joinOverlappingComponents = true;
          if (minAllowedDistance < 0)
              JMesh::error("minAllowedDistance must be >= 0.\n");
          else
              i++;
      }
  }*/
  else if (!strcmp(argv[i], "-u")) {
      if (i>=argc-1 || (uniformRemeshSteps = atoi(argv[i+1]))<1)
          JMesh::error("# uniform remesh steps must be >= 1.\n");
      i++;
  }
  else if (!strcmp(argv[i], "--decouple")) haveDecouple = true;
  else if (!strcmp(argv[i], "--no-clean")) clean = false;
  else if (!strcmp(argv[i], "-jc")) haveJoinClosestComponents = true;
  else if (!strcmp(argv[i], "-o")) {
      if (i<argc-1) {
          haveOutputFile = true;
          outputFile = argv[i+1];
          i++;
      }
  }
  else if (argv[i][0] == '-') JMesh::warning("%s - Unknown operation.\n",argv[i]);
 }

 // The loader performs the conversion to a set of oriented manifolds
 if (tin.load(argv[1]) != 0) JMesh::error("Can't open file.\n");
 // Join the second input argument if existing
 if (tin.append(argv[2]) == 0)
     JMesh::info("Joining the two meshfiles %s %s.\n", argv[1], argv[2]);
 input_filename = argv[1];

 // Keep only the biggest components
 tin.removeSmallestComponents( numberComponentsToKeep );

 // Fill holes by taking into account both sampling density and normal field continuity
 tin.fillSmallBoundaries(tin.E.numels(), true, true);

 if (joinOverlappingComponents) {
     tin.removeSmallestComponents(2);
     tin.joinOverlappingComponentPair();
//     tin.joinCloseOrOverlappingComponents( minAllowedDistance );
 }
 if (haveJoinClosestComponents)
 {
  printf("\nJoining input components ...\n");
  JMesh::begin_progress();
  while (joinClosestComponents(&tin, false, true, true)) JMesh::report_progress("Num. components: %d       ",tin.shells());
  JMesh::end_progress();
  tin.deselectTriangles();
 }

 if (uniformRemeshSteps) {
     printf("Uniform remeshing ...\n");
     tin.uniformRemesh(uniformRemeshSteps, 0, tin.E.numels());
 }

 if (haveDecouple) {
     printf("Decoupling first component from second one.\n");
     if(tin.shells() != 2) JMesh::warning("Incorrect number von components, won't decouple. Having %d and should have 2.", tin.shells());
     else tin.decoupleSecondFromFirstComponent(1, 50);
 }

 // Run geometry correction
 if (clean) {
     if (tin.boundaries() || !tin.clean(10, 3, numberComponentsToKeep)) {
      fprintf(stderr,"MeshFix failed!\n");
      fprintf(stderr,"Please try manually using ReMESH v1.2 or later (http://remesh.sourceforge.net).\n");
      FILE *fp = fopen("meshfix_log.txt","a");
      fprintf(fp,"MeshFix failed on %s\n",input_filename);
      fclose(fp);
     }
 }
 char *fname = createFilename( haveOutputFile ? outputFile : argv[1], subext, (save_vrml? ".wrl" : (save_stl? ".stl":".off")));
 printf("Saving output mesh to '%s'\n",fname);
 if (save_vrml)
     tin.saveVRML1(fname);
 else if (save_stl)
     tin.saveSTL(fname);
 else
     tin.saveOFF(fname);
 return 0;
}
