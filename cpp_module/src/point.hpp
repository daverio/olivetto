#ifndef POINT_HPP
#define POINT_HPP

#include <iostream>

struct Point{
  double x;
  double y;

  double norm();
};

Point operator+(const Point &a, const Point &b);
Point operator-(const Point &a, const Point &b);

Point operator*(const Point &p, const double &x);
Point operator/(const Point &p, const double &x);

Point& operator *=(Point& a, const double& b);
Point& operator /=(Point& a, const double& b);

double dist(const Point& a, const Point& b);

std::ostream& operator<<(std::ostream& os, const Point& dt);


#endif