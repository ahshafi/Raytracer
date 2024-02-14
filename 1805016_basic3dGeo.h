#include <math.h>
#ifndef BASIC3DGEO_H
#define BASIC3DGEO_H
const double EPS = 1e-5;
const double PI = acos(-1);
int dcmp(double x) { return abs(x) < EPS ? 0 : (x<0 ? -1 : 1);}
double degreeToRadian(double rad) { return rad*PI/180; }
struct Point {
    double x, y, z; Point() : x(0), y(0), z(0) {}
    Point(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    Point operator + (const Point& u) const {
        return Point(x + u.x, y + u.y, z + u.z); }
    Point operator - (const Point& u) const {
        return Point(x - u.x, y - u.y, z - u.z); }
    Point operator * (const double u) const {
        return Point(x * u, y * u, z * u); }
    Point operator / (const double u) const {
        return Point(x / u, y / u, z / u); }
    friend ostream &operator << (ostream &os,
         const Point &p) {return os<<p.x<<" "<<p.y << " " << p.z; }
};
double dot(Point a, Point b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
Point cross(Point a, Point b) { return 
    Point(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}
double length(Point a) { return sqrt(dot(a, a)); }
Point normalize(Point a) { return a/length(a); }
double distance(Point a, Point b) { return length(a-b); }
Point unit(const Point &p) { return p/length(p); }
// Rotate p around axis x, with angle radians.
Point rotate(Point p, Point axis, double angle) {
    axis = unit(axis); Point comp1 = p * cos(angle);
    Point comp2 = axis * (1-cos(angle)) * dot(axis, p);
    Point comp3 = cross(axis, p) * sin(angle);
    return comp1 + comp2 + comp3;
}

struct Line 
{
    Point a, v;
    Line(Point a, Point v): a(a), v(v) {}
};  ///a+tv
// returns the distance from point a to line l
double distancePointLine(Point p, Line l) {
    return length(cross(l.v, p - l.a)) / length(l.v); }
/// distance from Line ab to Line cd
double distanceLineLine(Line a, Line b) {
    Point cr = cross(a.v, b.v); double crl = length(cr);
    if (dcmp(crl) == 0)  return distancePointLine(a.a, b);
    return abs(dot(cr, a.a-b.a))/crl; }

struct Plane {
    Point normal;  /// Normal = (A, B, C)
    double d;      /// dot(Normal) = d <--> Ax + By + Cz = d
    Point P;      /// anyPoint on the plane, optional
    Plane(Point normal, double d) {
        double len = length(normal); assert(dcmp(len) > 0);
        normal = normal / len; d = d / len;
        if  (dcmp(normal.x))     P = Point(d/normal.x, 0, 0);
        else if (dcmp(normal.y)) P = Point(0, d/normal.y, 0);
        else                     P = Point(0, 0, d/normal.z);
    }
    ///Plane given by three Non-Collinear Points
    Plane(Point a, Point b, Point c)   {
        normal = unit(cross(b-a, c-a)); d = dot(normal,a); P=a;
    }
    bool onPlane(Point a) { return dcmp(dot(normal,a)-d)== 0;}
    double distance(Point a) { return abs(dot(normal, a) - d);}
    double isParallel(Line l){return dcmp(dot(l.v,normal))==0;}
//return t st l.a + t*l.v is a point on the plane, check 
//parallel first
    double intersectLine(Line l) {
        return dot(P-l.a, normal)/dot(l.v, normal); 
    } 
};
#endif