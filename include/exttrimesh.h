/****************************************************************************
* JMeshExt                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
*                                                                           *
* Copyright(C) 2006: IMATI-GE / CNR                                         *
*                                                                           *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef _EXTTRIMESH_H
#define _EXTTRIMESH_H

#include "jmesh.h"
#include "binTree.h"
#include "triangleoctree.h"
#include "component.h"
#include <set>
using std::set;

class ExtTriMesh : public Triangulation
{
 public:
 TriangleOctree *ot;
 // Constructors
 ExtTriMesh() : Triangulation() { ot = NULL; }
 ExtTriMesh(const char *s) : Triangulation(s) { ot = NULL; }
 ExtTriMesh(const Triangulation *t) : Triangulation(t) { ot = NULL; }
 ExtTriMesh(const Triangle *t, const bool keep_ref =false) : Triangulation(t, keep_ref) { ot = NULL; }

 Edge	*joinBoundaryLoops(bool =0, bool =1, bool =1); // (in "ALGORITHMS/holeFilling.C")
 Edge	*joinBoundaryLoops(Vertex *, Vertex *, bool =0, bool =1, bool =1); // (in "ALGORITHMS/holeFilling.C")
 int     fillSmallBoundaries(int, bool =0, bool =0);   // (in "ALGORITHMS/holeFilling.C")
 int     TriangulateHole(Edge *);		       // (in "ALGORITHMS/holeFilling.C")
 void    FillHole(Edge *, bool =0);		       // (in "ALGORITHMS/holeFilling.C")
 int refineSelectedHolePatches(Triangle * =NULL);      // (in "ALGORITHMS/holeFilling.C")
 void fairSelection(Triangle * =NULL);		       // (in "ALGORITHMS/holeFilling.C")

 // Mirko's functions
 void initializeOctree();
 int  joinCloseOrOverlappingComponents( double minAllowedDistance = 1.0 );
 int  joinOverlappingComponentPair();
 // returns the number of joined boundaries
 int joinComponentsCloseBoundaries(List *nl, List *ml, double maxDistanceToJoin);
 bool loopsHaveAllVerticesCloserThanDistance(List *loop, List *loop2, const double &distance);
 //! Determines the closest pair of vertices of the lists l1,l2. Good for boundary loops.
 double closestPair(List *l1, List *l2, Vertex **closest1, Vertex **closest2);
 //! Determines the closest partner of a vertex in a list of vertices
 double getClosestPartner(Vertex *v, List *l, Vertex **closestParnter);
 double mostDistantPartner(Vertex *v, List *l, Vertex **distantPartner);
 //! Returns true, if the Point p is inside the component. The component must be a
 //! closed surface. Searches for the closest vertex and uses the orientation of its triangles normal.
 //! Warning: The normals must be directed outwards of each component!
 bool isPointInComponent(Vertex *v, ComponentStruct *c);

 // Misc Algorithms (Implemented in "ALGORITHMS/*.C")

 void loopSubdivision(int);
 void modbutSubdivision();
 void sqrt3Subdivision();
 int laplacianSmooth(int =1, double =1.0);
 int uniformRemesh(int num_steps, int num_vertices = 0, int max_swaps = 10);
 int spherize(int);
 int featureRecover(double, double);
 int simplify(int, int =0, int =0, int =0);
 int multiplechoice_simplify(int, int =0, int =8, int =0);
 void mc_resample(int, int =0, int =0);
 int epsilonSample(double, int =0);

 int  selectIntersectingTriangles(UINT16 tri_per_cell=100);

 void tagPlanarRegionsBoundaries(double max_distance);

 //! Normalize all the shells and distribute them on a virtual sphere. O(N).
 void placeShellsOnVirtualSphere();
 Edge *flatten(Edge * =NULL);

 bool plumberSelect(Vertex *start, double radius);
};

#endif // _EXTTRIMESH_H
