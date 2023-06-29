#include "polygone_tests.hpp"


bool InsidePolygon(Point const * const polygon,size_t const N,Point const p)
{
  int counter = 0;
  double xinters;
  Point p1,p2;

  p1 = polygon[0];
  for (size_t i=1;i<=N;i++) {
    p2 = polygon[i % N];
    if (p.y > std::min(p1.y,p2.y)) {
      if (p.y <= std::max(p1.y,p2.y)) {
        if (p.x <= std::max(p1.x,p2.x)) {
          if (p1.y != p2.y) {
            xinters = (p.y-p1.y)*(p2.x-p1.x)/(p2.y-p1.y)+p1.x;
            if (p1.x == p2.x || p.x <= xinters)
              counter++;
          }
        }
      }
    }
    p1 = p2;
  }

  if (counter % 2 == 0)
    return false;
  else
    return true;
}

int orientation(const Point & A, const Point & B, const Point & C)
{
  double val = (B.y - A.y) * (C.x - B.x) - (B.x - A.x) * (C.y - B.y);
  if(val == 0) return 0;  // collinear
  return (val > 0)? 1: -1;     
}


bool check_colision(const Point & p1, const Point & q1, const Point & p2, const Point & q2,  Point & col)
{

    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    if (o1 != o2 && o3 != o4)
    {
      double A1 = q1.y - p1.y;
      double B1 = p1.x - q1.x;
      double C1 = A1 * p1.x + B1 * p1.y;
      double A2 = q2.y - p2.y;
      double B2 = p2.x - q2.x;
      double C2 = A2 * p2.x + B2 * p2.y;

      double det = A1 * B2 - A2 * B1;
      if(det == 0)
      {
        //colinear (parallel)
      }
      else
      {
        col.x = (B2 * C1 - B1 * C2) / det;
        col.y = (A1 * C2 - A2 * C1) / det;
        return true;
      }
    }

  return false;
}