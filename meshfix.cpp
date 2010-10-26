#include "exttrimesh.h"
#include <string.h>
#include <stdlib.h>

const char *input_filename;

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
 printf(" -ns <n>             Only the <n> biggest shells are kept.\n");
 printf(" -j                  Join 2 biggest components if they overlap, remove the rest.\n");
// printf(" -j <d>              Join components closer than <d> or overlapping.\n");
 printf(" -jc                 Join the closest pair of components.\n");
 printf(" -u <steps>          Uniform remeshing of the whole mesh, steps > 0\n");
 printf("   -nv <n>           Constrain number of vertices to <n> (only with -u)\n");
 printf(" --decouple <dmin>   Treat 1st file as inner, 2nd file as outer component.\n");
 printf("                     Resolve intersections by moving outer triangles outward.\n");
 printf("                     Constrain the distance between the components > dmin.\n");
 printf(" -d <dmin>           Synonym for --decouple <dmin>.\n");
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
 int uniformRemeshSteps = 0, numberOfVertices = 0;
 double decoupleMinDist = -1;
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
  else if (!strcmp(argv[i], "-nv")) {
      if (i>=argc-1 || (numberOfVertices = atoi(argv[i+1]))<1)
          JMesh::error("# of vertices must be >= 0.\n");
      i++;
  }
  else if (!strcmp(argv[i], "--decouple") || !strcmp(argv[i], "-d")) {
      if (i<argc-1) {
          decoupleMinDist = atof(argv[i+1]);
          if (decoupleMinDist < 0)
              JMesh::error("decoupleMinDist must be >= 0.\n");
          else
              i++;
      }
  }
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
     tin.uniformRemesh(uniformRemeshSteps, numberOfVertices, tin.E.numels());
 } else if(numberOfVertices) { JMesh::warning("-nv works only together with -u."); }

 if (decoupleMinDist >= 0) {
     printf("Decoupling first component from second one using %g as minimum allowed distance.\n", decoupleMinDist);
     if(tin.shells() != 2) JMesh::warning("Incorrect number of components, won't decouple. Having %d and should have 2.\n", tin.shells());
     else tin.decoupleSecondFromFirstComponent(decoupleMinDist, 20);
 }

 // Run geometry correction
 if (clean) {
     printf("Cleaning intersections, degeneracies ...\n");
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
