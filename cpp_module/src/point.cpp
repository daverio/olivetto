#include "point.hpp"
#include <cmath>

Point operator+(const Point &a, const Point &b)
{
    Point res = a;
    res.x += b.x;
    res.y += b.y;
    return res;
}   
Point operator-(const Point &a, const Point &b)
{
    Point res = a;
    res.x -= b.x;
    res.y -= b.y;
    return res;
}


Point operator*(const Point &p, const double &x)
{
    Point res = p;
    res.x *= x;
    res.y *= x;
    return res;
}
Point operator/(const Point &p, const double &x)
{
    Point res = p;
    res.x /= x;
    res.y /= x;
    return res;
}

Point& operator *=(Point& a, const double& b)
{
    a.x *= b;
    a.y *= b;
    return a;
}
Point& operator /=(Point& a, const double& b)
{
    a.x /= b;
    a.y /= b;
    return a;
}

Point& operator -=(Point& a, const Point& b)
{
    a.x -= b.x;
    a.y -= b.y;
    return a;
}
Point& operator +=(Point& a, const Point& b)
{
    a.x += b.x;
    a.y += b.y;
    return a;
}


double Point::norm()
{
    return std::sqrt(x*x+y*y);
}

double dist(const Point& a, const Point& b)
{
    return (b-a).norm();
}


std::ostream& operator<<(std::ostream& os, const Point& dt)
{
    os << '<' <<dt.x << ',' << dt.y << '>';
    return os;
}