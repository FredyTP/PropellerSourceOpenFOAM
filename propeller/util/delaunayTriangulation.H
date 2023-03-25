#ifndef PROPELLER_UTIL_DELAUNAY_TRIANGULATION
#define PROPELLER_UTIL_DELAUNAY_TRIANGULATION

#define EPSILON VSMALL

#include "List.H"
#include "label.H"
#include "point.H"
#include <functional>
namespace Foam
{

class delaunayTriangulation
{
public:
   typedef FixedList<label,3> ITRIANGLE;
   typedef FixedList<label,2> IEDGE;


/**
 * Construct voronoid diagram for all points inside the region 
 * with the function isInRegion
 * Cells are respect to each point provided, outside cells will be empty
 * Not all vertex are used by the cells
*/
static bool Voronoid
(
   const List<point>& points,
   List<point>& vertex, 
   List<List<label>> &cells,
   std::function<bool (point)> isInRegion,
   std::function<void (point,point,point&,bool)> intersectRegion
)
{
   List<point> nodeP = points;

   List<ITRIANGLE> triangles;
   delaunayTriangulation::Triangulate(nodeP,triangles,isInRegion,true);
   delaunayTriangulation::VoronoidFromTri(nodeP,triangles,vertex,cells,isInRegion,intersectRegion);
   return true;
}

static List<label> sortCounterClockwise(const List<point>& points, List<label>& poli,const point center)
{

   List<scalar> angles(poli.size());
   List<label> policopy = poli;
   for(label j = 0; j<poli.size();j++)
   {
         angles[j]=atan2(points[poli[j]].x() - center.x(),points[poli[j]].y()-center.y());
   }
   List<label> idx(poli.size());
   std::iota(idx.begin(), idx.end(), 0);
   //Sort index based on angles
   std::sort(idx.begin(),idx.end(),[angles](label a,label b)
   {
         return angles[a]<angles[b];
   });

   for(label j = 0; j<poli.size();j++)
   {
       poli[j]=policopy[idx[j]];
   }
   return idx;
}
static std::function<bool (point)> circularRegion(scalar radius)
{
   std::function<bool (point)> inCircle(
   [radius](point p1)
   {
      return magSqr(p1) <= sqr(radius);
   });
   return inCircle;
}
static std::function<void (point,point,point&,bool)> intersectCircle(scalar radius)
{

   std::function<void (point,point,point&,bool)> intersect(
      [radius](point p1, point p2, point& s1, bool p1_out)
   {
      scalar dx,dy,dr,D,sigdy,r;
      scalar xi1,xi2,yi1,yi2,xi,yi;

      scalar x1,x2,y1,y2;
      x1=p1.x();
      x2=p2.x();
      y1=p1.y();
      y2=p2.y();

      r=radius;

      dx=x2-x1;
      dy=y2-y1;
      dr=sqrt(pow(dx,2)+pow(dy,2));
      D=x1*y2-x2*y1;
      sigdy = dy<0?-1:1;
      xi1 = D*dy + sigdy*dx*sqrt(pow(r,2)*pow(dr,2)-pow(D,2));
      xi1/= pow(dr,2);

      xi2 = D*dy - sigdy*dx*sqrt(pow(r,2)*pow(dr,2)-pow(D,2));
      xi2/= pow(dr,2);

      yi1 = -D*dx + std::abs(dy)*sqrt(pow(r,2)*pow(dr,2)-pow(D,2));
      yi1/= pow(dr,2);

      yi2 = -D*dx - std::abs(dy)*sqrt(pow(r,2)*pow(dr,2)-pow(D,2));
      yi2/= pow(dr,2);

      scalar xout,yout;
      xout = p1_out?x1:x2;
      yout = p1_out?y1:y2;

      scalar d1 = pow(xout-xi1,2)+pow(yout-yi1,2);
      scalar d2 = pow(xout-xi2,2)+pow(yout-yi2,2);

      if(d1<d2)
      {
         xi=xi1;
         yi=yi1;
      }
      else
      {
         xi=xi2;
         yi=yi2;
      }

      s1.x()=xi;
      s1.y()=yi;
   });
   return intersect;
}


static int VoronoidFromTri
(
   List<point>& nucli,
   List<ITRIANGLE>& tri, 
   List<point>& orto, 
   List<List<label>>& voronoid,
   std::function<bool (point)> isInRegion,
   std::function<void (point,point,point&,bool)> findIntersection
)
{

   orto.resize(tri.size());
   /*Compute orto centers*/
   for(label i=0;i <tri.size();i++)
   {
      double xc,yc,r;
      double x1,y1,x2,y2,x3,y3;
      
      //Get triangle vertex
      x1 = nucli[tri[i][0]].x();
      y1 = nucli[tri[i][0]].y();

      x2 = nucli[tri[i][1]].x();
      y2 = nucli[tri[i][1]].y();

      x3 = nucli[tri[i][2]].x();
      y3 = nucli[tri[i][2]].y();

      delaunayTriangulation::CircumCircle(0,0,x1,y1,x2,y2,x3,y3,&xc,&yc,&r);

      //list of triangle ortocenters
      orto[i]= point{xc,yc,0};
   }

   //List of index for each nucli
   voronoid.resize(nucli.size());

   //Traverse all triangles and get index
   for(label i=0; i <tri.size();i++)
   {
      voronoid[tri[i][0]].append(i);
      voronoid[tri[i][1]].append(i);
      voronoid[tri[i][2]].append(i);
   }
  
   for(label i=0;i<voronoid.size();i++)
   {
      // Sort voronoid points in conterclockwise direction
      sortCounterClockwise(orto,voronoid[i],nucli[i]);

      //Create temp poligon
      List<label> poli = voronoid[i];
      //Add points inside domain or intersection with domain
      List<label> insidePoli;

      //List of visited points
      List<bool> visited(poli.size(),false);

      //If point is outside continue
      if(!isInRegion(nucli[i]))
      {
         voronoid[i].clear();
         continue;
      }
      for(label j = 0; j<poli.size();j++)
      {
            
            point po;
            po = orto[poli[j]];               
            if(!isInRegion(po))
            {
               if(!visited[j])
               {
                  point p1,p2,int1,int2,pout1,pout2;
                  label i1,i2,o1,o2;
                  
                  List<label> to_remove;

                  i1=j;
                  i2=j;
                  do
                  {
                     visited[i1]=true;
                     o1=i1;
                     i1 = i1==0?poli.size()-1 : i1-1;
                     p1 = orto[poli[i1]];
                  } while (!isInRegion(p1));
                  
                  //Get next point indexo1 = (i1+1)%poli.size();

                  do
                  {
                     visited[i2]=true;
                     o2=i2;
                     i2 = (i2+1)%poli.size();
                     p2 = orto[poli[i2]];

                  } while (!isInRegion(p2));
                  
                  
               
                  if(o1 == o2)
                  {
                     findIntersection(orto[poli[o1]],p1,int1,true);
                     findIntersection(orto[poli[o1]],p2,int2,true);
                     insidePoli.append(orto.size());
                     orto.append(int1);
                     insidePoli.append(orto.size());
                     orto.append(int2);
                  }
                  else
                  {
                     pout1 = orto[poli[o1]];
                     pout2 = orto[poli[o2]];

                     findIntersection(pout1,p1,int1,true);
                     findIntersection(pout2,p2,int2,true);

                     insidePoli.append(orto.size());
                     orto.append(int1);

                     insidePoli.append(orto.size());
                     orto.append(int2);
                  }
                  
               }
            }
            else
            {

               insidePoli.append(poli[j]);
               
            }
           
      }
      voronoid[i] = insidePoli;
   }

   //Last 3 points are from supertriangle
   nucli.resize(nucli.size()-3);
   voronoid.resize(nucli.size());

   return 0;
}

/**
 * Remove all triangles with index >= nv
*/
static void removeSuperTri(List<ITRIANGLE>& tri,label nv)
{
   label ntri = tri.size();
   for (label i=0;i<ntri;i++) 
   {
      if (tri[i][0] >= nv || tri[i][1] >= nv || tri[i][2] >= nv) 
      {
         tri[i] = tri[tri.size()-1];
         i--;
         ntri--;
      }
   }
   tri.resize(ntri);
}
/*
  Implemented from http://paulbourke.net/papers/triangulate/

   Triangulation subroutine
   Takes as input NV vertices in array ppoint
   Returned is a list of ntri triangular faces in the array v
   These triangles are arranged in a consistent clockwise order.
   The triangle array 'v' should be malloced to 3 * nv
   The vertex array ppoint must be big enough to hold 3 more points
   The vertex array must be sorted in increasing x values say

   qsort(p,nv,sizeof(point),pointCompare);
      :
   int pointCompare(void *v1,void *v2)
   {
      point *p1,*p2;
      p1 = v1;
      p2 = v2;
      if (p1->x < p2->x)
         return(-1);
      else if (p1->x > p2->x)
         return(1);
      else
         return(0);
   }
*/
static int Triangulate(List<point>& pxyz,List<ITRIANGLE>& v,std::function<bool (point)> isInRegion = [](point p){return false;},bool keepSuperTri = false)
{
   List<bool> complete;
   List<IEDGE> edges;
   int nedge = 0;
   int trimax,emax = 200;
   int status = 0;

   int inside;
   int i,j,k;
   double xp,yp,x1,y1,x2,y2,x3,y3,xc,yc,r;
   double xmin,xmax,ymin,ymax,xmid,ymid;
   double dx,dy,dmax;

   //Original number of vertex
   label nv = pxyz.size();
   label ntri = 0;

   //resize to allocate the 3 supertriangle vertex
   pxyz.resize(nv+3);

   /* Allocate memory for the completeness list, flag for each triangle */
   trimax = 4 * nv;

   //Set vector to check complete triangles
   complete.resize(trimax);
   v.resize(trimax);
   
   /*if ((complete = new int[trimax]) == NULL) {
      status = 1;
      goto skip;
   }*/

   /* Allocate memory for the edge list */
   edges.resize(emax);
   /*if ((edges = new IEDGE[emax]) == NULL) {
      status = 2;
      goto skip;
   }*/

   /*
      Find the maximum and minimum vertex bounds.
      This is to allow calculation of the bounding triangle
   */
   xmin = pxyz[0].x();
   ymin = pxyz[0].y();
   xmax = xmin;
   ymax = ymin;
   for (i=1;i<nv;i++) {
      if (pxyz[i].x() < xmin) xmin = pxyz[i].x();
      if (pxyz[i].x() > xmax) xmax = pxyz[i].x();
      if (pxyz[i].y() < ymin) ymin = pxyz[i].y();
      if (pxyz[i].y() > ymax) ymax = pxyz[i].y();
   }

   scalar factor=1;
   do
   {
      dx = xmax - xmin;
      dy = ymax - ymin;
      dmax = factor*((dx > dy) ? dx : dy);
      xmid = (xmax + xmin) / 2.0;
      ymid = (ymax + ymin) / 2.0;

      pxyz[nv+0].x() = xmid - 20 * dmax;
      pxyz[nv+0].y() = ymid - dmax;
      pxyz[nv+0].z() = 0.0;
      pxyz[nv+1].x() = xmid;
      pxyz[nv+1].y() = ymid + 20 * dmax;
      pxyz[nv+1].z() = 0.0;
      pxyz[nv+2].x() = xmid + 20 * dmax;
      pxyz[nv+2].y() = ymid - dmax;
      pxyz[nv+2].z() = 0.0;

      factor++;
   } while (isInRegion(pxyz[nv+0]) || isInRegion(pxyz[nv+1]) || isInRegion(pxyz[nv+2]));
   


   /*
      Set up the supertriangle
      This is a triangle which encompasses all the sample points.
      The supertriangle coordinates are added to the end of the
      vertex list. The supertriangle is the first triangle in
      the triangle list.
   */
  
   v[0][0] = nv;
   v[0][1] = nv+1;
   v[0][2] = nv+2;
   complete[0] = false;
   ntri = 1;

   /*
      Include each point one at a time into the existing mesh
   */
   for (i=0;i<nv;i++) {

      xp = pxyz[i].x();
      yp = pxyz[i].y();
      nedge = 0;

      /*
         Set up the edge buffer.
         If the point (xp,yp) lies inside the circumcircle then the
         three edges of that triangle are added to the edge buffer
         and that triangle is removed.
      */
      for (j=0;j<ntri;j++) {
         if (complete[j])
            continue;
         x1 = pxyz[v[j][0]].x();
         y1 = pxyz[v[j][0]].y();
         x2 = pxyz[v[j][1]].x();
         y2 = pxyz[v[j][1]].y();
         x3 = pxyz[v[j][2]].x();
         y3 = pxyz[v[j][2]].y();
         inside = CircumCircle(xp,yp,x1,y1,x2,y2,x3,y3,&xc,&yc,&r);
         if (xc < xp && ((xp-xc)*(xp-xc)) > r)
				complete[j] = true;
         if (inside) {
            /* Check that we haven't exceeded the edge list size */
            if (nedge+3 >= emax) {
               emax += 100;
               edges.resize(emax);
            }
            edges[nedge+0][0] = v[j][0];
            edges[nedge+0][1] = v[j][1];
            edges[nedge+1][0] = v[j][1];
            edges[nedge+1][1] = v[j][2];
            edges[nedge+2][0] = v[j][2];
            edges[nedge+2][1] = v[j][0];
            nedge += 3;
            v[j] = v[ntri-1];
            complete[j] = complete[ntri-1];
            ntri--;
            j--;
         }
      }

      /*
         Tag multiple edges
         Note: if all triangles are specified anticlockwise then all
               interior edges are opposite pointing in direction.
      */
      for (j=0;j<nedge-1;j++) {
         for (k=j+1;k<nedge;k++) {
            if ((edges[j][0] == edges[k][1]) && (edges[j][1] == edges[k][0])) {
               edges[j][0] = -1;
               edges[j][1] = -1;
               edges[k][0] = -1;
               edges[k][1] = -1;
            }
            /* Shouldn't need the following, see note above */
            if ((edges[j][0] == edges[k][0]) && (edges[j][1] == edges[k][1])) {
               edges[j][0] = -1;
               edges[j][1] = -1;
               edges[k][0] = -1;
               edges[k][1] = -1;
            }
         }
      }

      /*
         Form new triangles for the current point
         Skipping over any tagged edges.
         All edges are arranged in clockwise order.
      */
      for (j=0;j<nedge;j++) {
         if (edges[j][0] < 0 || edges[j][1] < 0)
            continue;
         if (ntri >= trimax) {
            status = 4;
            return(status);
         }
         v[ntri][0] = edges[j][0];
         v[ntri][1] = edges[j][1];
         v[ntri][2] = i;
         complete[ntri] = false;
         ntri++;
      }
   }

   v.resize(ntri);
   /*
      Remove triangles with supertriangle vertices
      These are triangles which have a vertex number greater than nv
   */
  if(!keepSuperTri)
  {
      delaunayTriangulation::removeSuperTri(v,nv);
      pxyz.resize(nv);
  }

  return(status);
   
}
/*
   Return true if a point (xp,yp) is inside the circumcircle made up
   of the points (x1,y1), (x2,y2), (x3,y3)
   The circumcircle centre is returned in (xc,yc) and the radius r
   NOTE: A point on the edge is inside the circumcircle
*/
static int CircumCircle(double xp,double yp,
   double x1,double y1,double x2,double y2,double x3,double y3,
   double *xc,double *yc,double *rsqr)
{
   double m1,m2,mx1,mx2,my1,my2;
   double dx,dy,drsqr;
   double fabsy1y2 = fabs(y1-y2);
   double fabsy2y3 = fabs(y2-y3);

   /* Check for coincident points */
   if (fabsy1y2 < EPSILON && fabsy2y3 < EPSILON)
       return(false);

   if (fabsy1y2 < EPSILON) {
      m2 = - (x3-x2) / (y3-y2);
      mx2 = (x2 + x3) / 2.0;
      my2 = (y2 + y3) / 2.0;
      *xc = (x2 + x1) / 2.0;
      *yc = m2 * (*xc - mx2) + my2;
   } else if (fabsy2y3 < EPSILON) {
      m1 = - (x2-x1) / (y2-y1);
      mx1 = (x1 + x2) / 2.0;
      my1 = (y1 + y2) / 2.0;
      *xc = (x3 + x2) / 2.0;
      *yc = m1 * (*xc - mx1) + my1;
   } else {
      m1 = - (x2-x1) / (y2-y1);
      m2 = - (x3-x2) / (y3-y2);
      mx1 = (x1 + x2) / 2.0;
      mx2 = (x2 + x3) / 2.0;
      my1 = (y1 + y2) / 2.0;
      my2 = (y2 + y3) / 2.0;
      *xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
      if (fabsy1y2 > fabsy2y3) {
         *yc = m1 * (*xc - mx1) + my1;
      } else {
         *yc = m2 * (*xc - mx2) + my2;
      }
   }

   dx = x2 - *xc;
   dy = y2 - *yc;
   *rsqr = dx*dx + dy*dy;

   dx = xp - *xc;
   dy = yp - *yc;
   drsqr = dx*dx + dy*dy;

   // Original
   //return((drsqr <= *rsqr) ? true : false);
   // Proposed by Chuck Morris
   return((drsqr - *rsqr) <= EPSILON ? true : false);
}


};


}

#endif /* PROPELLER_UTIL_DELAUNAY_TRIANGULATION */




