#ifndef POLIGONE_TESTS_HPP
#define POLIGONE_TESTS_HPP

#include "point.hpp"

bool InsidePolygon(Point const * const polygon,size_t const N,Point const p);

int orientation(const Point & A, const Point & B, const Point & C);

bool check_colision(const Point & p1, const Point & q1, const Point & p2, const Point & q2,  Point & col);



#endif