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