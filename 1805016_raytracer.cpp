#include<bits/stdc++.h>
using namespace std;
#include <GL/glut.h>
#define pi (2 * acos(0.0))
#include "1805016_basic3dGeo.h"
#include "1805016_debug.h"
#include "1805016_bitmap_image.hpp"
#define dbg1(args...) //printf(args)
double cameraHeight;
double cameraAngle;
double cameraRadius;
Point cameraPos;
Point cameraRef;
Point cameraUp;
double angle;

double near, far;
double fovY;
double aspectRatio;
int recursionLevel;
int resolution;
int objectCount;
//vector<Pyramid> pyramids;
//vector<Cube> cubes;

class CheckerBoard;
class Spherre;
class Cube;
class NormalLight;

void generateRays();
bool drawgrid;
void drawGrid()
{
    int i;
    if (drawgrid == 1)
    {
        glColor3f(0, 1, 0); // grey
        glBegin(GL_LINES);
        glVertex3f(-100, 0, 10);
        glVertex3f(100, 0, 10);
        glVertex3f(0, -100, 10);
        glVertex3f(0, 100, 10);
        glVertex3f(0, 0, -100);
        glVertex3f(0, 0, 100);

        glColor3f(.5, .5, 0); // grey
        for (i = -100; i <= 100; i++)
        {

            if (i == 0)
                continue; // SKIP the MAIN axes

            // lines parallel to Y-axis
            glVertex3f(i * 10, -1000, 10);
            glVertex3f(i * 10, 1000, 10);

            // lines parallel to X-axis
            glVertex3f(-1000, i * 10, 10);
            glVertex3f(1000, i * 10, 10);
        }

        glEnd();
    }
}
void keyboardListener(unsigned char key, int x, int y)
{
    const double SPEED = .3;
    Point cameraSideDir = normalize(cross(cameraRef - cameraPos, cameraUp));
    Point cameraForwardDir = normalize(cameraRef - cameraPos);
    switch (key)
    {
        case '0':
            generateRays();
            break;
        case '3':
            {cameraUp = normalize(rotate(cameraUp, cameraSideDir, SPEED));
            cameraForwardDir = rotate(cameraForwardDir, cameraSideDir, SPEED);
            cameraRef = cameraPos + cameraForwardDir;
            break;}
        case '4':
            {cameraUp = normalize(rotate(cameraUp, cameraSideDir, -SPEED));
            cameraForwardDir = rotate(cameraForwardDir, cameraSideDir, -SPEED);
            cameraRef = cameraPos + cameraForwardDir;
            break;}
        case '1':
            {cameraForwardDir = normalize(rotate(cameraForwardDir, cameraUp, SPEED));
            cameraRef = cameraPos + cameraForwardDir;
            break;}
        case '2':
            {cameraForwardDir = normalize(rotate(cameraForwardDir, cameraUp, -SPEED));
            cameraRef = cameraPos + cameraForwardDir;
            break;}
        case '5':
            cameraUp = normalize(rotate(cameraUp, cameraForwardDir, -SPEED));
            break;
        case '6':
            cameraUp = normalize(rotate(cameraUp, cameraForwardDir, SPEED));
            break;
        case ',':
            drawgrid = !drawgrid;
            break;
        case '.':
            break;
        case 'a':
            angle -= 2;
            break;
        case 'd':
            angle += 2;
            break;
        default:
            break;
    }
}

void specialKeyListener(int key, int x, int y)
{
    Point cameraForwardDir = (cameraRef - cameraPos) / length(cameraRef - cameraPos);
    Point cameraUpDir = cameraUp / length(cameraUp);
    Point cameraSideDir = cross(cameraForwardDir, cameraUpDir) / length(cross(cameraForwardDir, cameraUpDir));
    const double SPEED = 10;
    switch (key)
    {
        case GLUT_KEY_DOWN: // down arrow key
            cameraPos = cameraPos - cameraForwardDir * SPEED;
            cameraRef = cameraRef - cameraForwardDir * SPEED;
            break;
        case GLUT_KEY_UP: // up arrow key
            cameraPos = cameraPos + cameraForwardDir * SPEED;
            cameraRef = cameraRef + cameraForwardDir * SPEED;
            break;

        case GLUT_KEY_RIGHT:
            cameraPos = cameraPos + cameraSideDir * SPEED;
            cameraRef = cameraRef + cameraSideDir * SPEED;
            break;
        case GLUT_KEY_LEFT:
            cameraPos = cameraPos - cameraSideDir * SPEED;
            cameraRef = cameraRef - cameraSideDir * SPEED;
            break;

        case GLUT_KEY_PAGE_UP:
            cameraPos = cameraPos + cameraUpDir * SPEED;
            cameraRef = cameraRef + cameraUpDir * SPEED;
            break;
        case GLUT_KEY_PAGE_DOWN:
            cameraPos = cameraPos - cameraUpDir * SPEED;
            cameraRef = cameraRef - cameraUpDir * SPEED;
            break;

        case GLUT_KEY_INSERT:
            break;

        case GLUT_KEY_HOME:
            break;
        case GLUT_KEY_END:
            break;

        default:
            break;
    }
}

void mouseListener(int button, int state, int x, int y)
{ // x, y is the x-y of the screen (2D)
    switch (button)
    {
    case GLUT_LEFT_BUTTON:
        if (state == GLUT_DOWN)
        { // 2 times?? in ONE click? -- solution is checking DOWN or UP
            //drawaxes = 1 - drawaxes;
        }
        break;

    case GLUT_RIGHT_BUTTON:
        //........
        break;

    case GLUT_MIDDLE_BUTTON:
        //........
        break;

    default:
        break;
    }
}

int boardWidth = 100; // Number of cells in a row/column

class Ray
{
public:
    Point start;
    Point dir;
    string object;
    Ray(Point start, Point dir) : start(start), dir(dir) {}
};

class CheckerBoard
{
public:
    int cellSize; // Size of each cell in pixels
    Point normal = Point(0, 0, 1);
    Point center = Point(0, 0, 0);
    double ambient, diffuse, reflection;
    void drawCheckerboard(int boardWidth) {

        for (int i = -boardWidth + cameraPos.x / cellSize; i < boardWidth + cameraPos.x / cellSize; ++i) {
            for (int j = -boardWidth + cameraPos.y / cellSize; j < boardWidth + cameraPos.y / cellSize; ++j) {
                if ((i + j) % 2 == 0) {
                    glColor3f(1.0f, 1.0f, 1.0f); // White cell
                } else {
                    glColor3f(0.0f, 0.0f, 0.0f); // Black cell
                }

                glBegin(GL_QUADS);
                glVertex2i(i * cellSize, j * cellSize);
                glVertex2i((i + 1) * cellSize, j * cellSize);
                glVertex2i((i + 1) * cellSize, (j + 1) * cellSize);
                glVertex2i(i * cellSize, (j + 1) * cellSize);
                glEnd();
            }
        }

    }
    bool doesIntersect(Ray ray)
    {
        double t = dot(normal, center - ray.start) / dot(normal, ray.dir);
        if(t < 0) return false;
        else return true;
    }
    Point intersection(Ray ray)
    {
        double t = dot(normal, center - ray.start) / dot(normal, ray.dir);
        return ray.start + ray.dir * t;
    }
    Point getColor(Point p)
    {
        int i = floor((p.x - center.x) / cellSize);
        int j = floor((p.y - center.y) / cellSize);
        //cout << i << " " << j << endl;
        if((i + j) % 2 == 0)
        {
            return Point(1, 1, 1);
        }
        else
        {
            return Point(0, 0, 0);
        }
    }
}checkerBoard;

class Spherre
{
public:
    Point center;
    double radius;
    double r, g, b;
    double ambient, diffuse, specular, reflection;
    double shininess;
    Spherre() {};
    Spherre(Point center, double radius, double r, double g, double b, double ambient, double diffuse, double specular, double reflection, double shininess) : center(center), radius(radius), r(r), g(g), b(b), ambient(ambient), diffuse(diffuse), specular(specular), reflection(reflection), shininess(shininess) {}

    void drawSphereSegments()
    {
        glColor3f(r, g, b);
        glPushMatrix();
        glTranslatef(center.x, center.y, center.z);
        glutSolidSphere(radius, 100, 100);
        glPopMatrix();
    }
    bool insideSphere(Point p)
    {
        return distance(p, center) <= radius;
    }
    bool doesIntersect(Ray ray)
    {
        if(insideSphere(ray.start)) return false;

        Point L = center - ray.start;
        double tca = dot(L, ray.dir);
        if(tca < 0) return false;
        double d2 = dot(L, L) - tca * tca;
        if(d2 > radius * radius) return false;
        
        return true;
    }

    Point pointOfIntersection(Ray ray)
    {
        Point L = center - ray.start;
        double tca = dot(L, ray.dir); // tca = t closest approach
        double d2 = dot(L, L) - tca * tca;
        double thc = sqrt(radius * radius - d2); // thc = t hit circle
        double d = tca - thc;
        return ray.start + ray.dir * d;
    }

    Point normalAtPoint(Point point)
    {
        return normalize(point - center);
    }
};
vector<Spherre> spheres;

class Pyramid
{
public:
    Point bottomLeft;
    double side;
    double height;
    double r, g, b;
    double ambient, diffuse, specular, reflection;
    double shininess;
    Point curNormal;
    vector<array<Point, 3>> triangles;
    Point base[4];
    Pyramid() {};
    Pyramid(Point bottomLeft, double side, double height, double r, double g, double b, double ambient, double diffuse, double specular, double reflection, double shininess) : bottomLeft(bottomLeft), side(side), height(height), r(r), g(g), b(b), ambient(ambient), diffuse(diffuse), specular(specular), reflection(reflection), shininess(shininess) 
    {
        double centerX = bottomLeft.x + side / 2, centerY = bottomLeft.y + side / 2, centerZ = bottomLeft.z;
        Point center(centerX, centerY, centerZ);
        base[0] = center + Point(-side / 2, -side / 2, 0);
        base[1] = center + Point(side / 2, -side / 2, 0);
        base[2] = center + Point(side / 2, side / 2, 0);
        base[3] = center + Point(-side / 2, side / 2, 0);
        triangles.push_back({center + Point(-side / 2, -side / 2, 0), center + Point(side / 2, -side / 2, 0), center + Point(0, 0, height)});
        triangles.push_back({center + Point(side / 2, -side / 2, 0), center + Point(side / 2, side / 2, 0), center + Point(0, 0, height)});
        triangles.push_back({center + Point(side / 2, side / 2, 0), center + Point(-side / 2, side / 2, 0), center + Point(0, 0, height)});
        triangles.push_back({center + Point(-side / 2, side / 2, 0), center + Point(-side / 2, -side / 2, 0), center + Point(0, 0, height)});
    }

    void drawPyramid()
    {
        glColor3f(r, g, b);
        glPushMatrix();
        glTranslatef(bottomLeft.x + side / 2, bottomLeft.y + side / 2, bottomLeft.z);
        glBegin(GL_QUADS);
        {
            glVertex3f(-side / 2, -side / 2, 0);
            glVertex3f(side / 2, -side / 2, 0);
            glVertex3f(side / 2, side / 2, 0);
            glVertex3f(-side / 2, side / 2, 0);
        }
        glEnd();
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(-side / 2, -side / 2, 0);
            glVertex3f(side / 2, -side / 2, 0);
            glVertex3f(0, 0, height);
        }
        glEnd();
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(side / 2, -side / 2, 0);
            glVertex3f(side / 2, side / 2, 0);
            glVertex3f(0, 0, height);
        }
        glEnd();
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(side / 2, side / 2, 0);
            glVertex3f(-side / 2, side / 2, 0);
            glVertex3f(0, 0, height);
        }
        glEnd();
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(-side / 2, side / 2, 0);
            glVertex3f(-side / 2, -side / 2, 0);
            glVertex3f(0, 0, height);
        }
        glEnd();
        glPopMatrix();
    }

    bool insidePyramid(Point p)
    {
        Point center = bottomLeft + Point(side / 2, side / 2, 0);
        Point normal = Point(0, 0, 1);
        Point A = center + Point(-side / 2, -side / 2, 0);
        Point B = center + Point(side / 2, -side / 2, 0);
        Point C = center + Point(side / 2, side / 2, 0);
        Point D = center + Point(-side / 2, side / 2, 0);
        Point AB = B - A;
        Point AC = C - A;
        Point AD = D - A;
        Point AP = p - A;
        Point BC = C - B;
        Point BP = p - B;
        Point CD = D - C;
        Point CP = p - C;
        Point DP = p - D;
        Point normal1 = cross(AB, AP);
        Point normal2 = cross(BC, BP);
        Point normal3 = cross(CD, CP);
        Point normal4 = cross(DP, DP);
        if(dot(normal, normal1) >= 0 && dot(normal, normal2) >= 0 && dot(normal, normal3) >= 0 && dot(normal, normal4) >= 0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    bool doesIntersectBase(Ray ray)
    {
        if(insidePyramid(ray.start)) return false;
        for(int i = 0; i < 4; i++)
        {
            Point A = base[i];
            Point B = base[(i + 1) % 4];
            Point C = base[(i + 2) % 4];
            Point AB = B - A;
            Point AC = C - A;
            Point normal = cross(AB, AC);
            double t = dot(normal, A - ray.start) / dot(normal, ray.dir);
            if(t < 0) continue;
            Point P = ray.start + ray.dir * t;
            Point AP = P - A;
            Point BP = P - B;
            Point CP = P - C;
            Point normal1 = cross(AB, AP);
            Point normal2 = cross(AC, AP);
            Point normal3 = cross(AC, CP);
            if(dot(normal, normal1) >= 0 && dot(normal, normal2) >= 0 && dot(normal, normal3) >= 0)
            {
                return true;
            }
        }

        return false;
    }

    bool doesIntersectTriangularFaces(Ray ray)
    {
        if(insidePyramid(ray.start)) return false;
        for(int i = 0; i < 4; i++)
        {
            Point A = triangles[i][0];
            Point B = triangles[i][1];
            Point C = triangles[i][2];
            Point AB = B - A;
            Point AC = C - A;
            Point normal = cross(AB, AC);
            double t = dot(normal, A - ray.start) / dot(normal, ray.dir);
            if(t < 0) continue;
            Point P = ray.start + ray.dir * t;
            Point AP = P - A;
            Point BP = P - B;
            Point CP = P - C;
            Point normal1 = cross(AB, AP);
            Point normal2 = cross(AC, AP);
            Point normal3 = cross(AC, CP);
            if(dot(normal, normal1) >= 0 && dot(normal, normal2) >= 0 && dot(normal, normal3) >= 0)
            {
                return true;
            }
        }

        return false;
    }

    bool doesIntersect(Ray ray)
    {
        return doesIntersectBase(ray) || doesIntersectTriangularFaces(ray);
    }

    Point pointOfIntersection(Ray ray)
    {
        double t = 1e9;
        Point point;
        if(doesIntersectBase(ray))
        {
            for(int i = 0; i < 4; i++)
            {
                Point A = base[i];
                Point B = base[(i + 1) % 4];
                Point C = base[(i + 2) % 4];
                Point AB = B - A;
                Point AC = C - A;
                Point normal = cross(AB, AC);
                double t1 = dot(normal, A - ray.start) / dot(normal, ray.dir);
                if(t1 < 0) continue;
                Point P = ray.start + ray.dir * t1;
                Point AP = P - A;
                Point BP = P - B;
                Point CP = P - C;
                Point normal1 = cross(AB, AP);
                Point normal2 = cross(AC, AP);
                Point normal3 = cross(AC, CP);
                if(dot(normal, normal1) >= 0 && dot(normal, normal2) >= 0 && dot(normal, normal3) >= 0)
                {
                    if(t1 < t)
                    {
                        t = t1;
                        point = P;
                        curNormal = Point(0, 0, 1);
                    }
                }
            }
        }
        if(doesIntersectTriangularFaces(ray))
        {
            for(int i = 0; i < 4; i++)
            {
                Point A = triangles[i][0];
                Point B = triangles[i][1];
                Point C = triangles[i][2];
                Point AB = B - A;
                Point AC = C - A;
                Point normal = cross(AB, AC);
                double t1 = dot(normal, A - ray.start) / dot(normal, ray.dir);
                if(t1 < 0) continue;
                Point P = ray.start + ray.dir * t1;
                Point AP = P - A;
                Point BP = P - B;
                Point CP = P - C;
                Point normal1 = cross(AB, AP);
                Point normal2 = cross(AC, AP);
                Point normal3 = cross(AC, CP);
                if(dot(normal, normal1) >= 0 && dot(normal, normal2) >= 0 && dot(normal, normal3) >= 0)
                {
                    if(t1 < t)
                    {
                        t = t1;
                        point = P;
                        curNormal = normalize(cross(A - C, B - C));

                    }
                }
            }
        }
        return point;
    }

    Point normalAtPonit(Point point)
    {
        return curNormal;
    }
};
vector<Pyramid> pyramids;

class Cube
{
public:
    Point bottomLeft;
    double side;
    double r, g, b;
    double ambient, diffuse, specular, reflection;
    double shininess;
    Cube() {};
    Cube(Point bottomLeft, double side, double r, double g, double b, double ambient, double diffuse, double specular, double reflection, double shininess) : bottomLeft(bottomLeft), side(side), r(r), g(g), b(b), ambient(ambient), diffuse(diffuse), specular(specular), reflection(reflection), shininess(shininess) {}

    void DrawCube()
    {
        glColor3f(r, g, b);
        glPushMatrix();
        glTranslatef(bottomLeft.x + side / 2, bottomLeft.y + side / 2, bottomLeft.z + side / 2);
        glutSolidCube(side);
        glPopMatrix();
    }

    Point normalAtPoint(Point point)
    {
        Point normal;
        if(dcmp(point.x - bottomLeft.x) == 0)
        {
            normal = Point(-1, 0, 0);
        }
        else if(dcmp(point.x - (bottomLeft.x + side)) == 0)
        {
            normal = Point(1, 0, 0);
        }
        else if(dcmp(point.y - bottomLeft.y))
        {
            normal = Point(0, -1, 0);
        }
        else if(dcmp(point.y - (bottomLeft.y + side)) == 0)
        {
            normal = Point(0, 1, 0);
        }
        else if(dcmp(point.z - bottomLeft.z) == 0)
        {
            normal = Point(0, 0, -1);
        }
        else if(dcmp(point.z - (bottomLeft.z + side)) == 0)
        {
            normal = Point(0, 0, 1);
        }
        return normal;
    }

    bool insideCube(Point p)
    {
        return p.x >= bottomLeft.x && p.x <= bottomLeft.x + side && p.y >= bottomLeft.y && p.y <= bottomLeft.y + side && p.z >= bottomLeft.z && p.z <= bottomLeft.z + side;
    }

    // bool doesIntersect(Ray ray)
    // {
    //     if(insideCube(ray.start)) return false;

    //     double t1 = (bottomLeft.x - ray.start.x) / ray.dir.x;
    //     double t2 = (bottomLeft.x + side - ray.start.x) / ray.dir.x;
    //     double t3 = (bottomLeft.y - ray.start.y) / ray.dir.y;
    //     double t4 = (bottomLeft.y + side - ray.start.y) / ray.dir.y;
    //     double t5 = (bottomLeft.z - ray.start.z) / ray.dir.z;
    //     double t6 = (bottomLeft.z + side - ray.start.z) / ray.dir.z;

    //     double tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
    //     double tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));

    //     if(tmax < 0) return false;
    //     if(tmin > tmax) return false;

    //     return true;
    // }

    // Point pointOfIntersection(Ray ray)
    // {
    //     double t1 = (bottomLeft.x - ray.start.x) / ray.dir.x;
    //     double t2 = (bottomLeft.x + side - ray.start.x) / ray.dir.x;
    //     double t3 = (bottomLeft.y - ray.start.y) / ray.dir.y;
    //     double t4 = (bottomLeft.y + side - ray.start.y) / ray.dir.y;
    //     double t5 = (bottomLeft.z - ray.start.z) / ray.dir.z;
    //     double t6 = (bottomLeft.z + side - ray.start.z) / ray.dir.z;

    //     double tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
    //     double tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));

    //     double t = tmin;

    //     return ray.start + ray.dir * t;
    // }

    bool doesIntersect(Ray ray)
    {
        Point center = bottomLeft + Point(side / 2, side / 2, side / 2);
        double halfSize = side / 2;
        Point minBounds = center - Point(halfSize, halfSize, halfSize);
        Point maxBounds = center + Point(halfSize, halfSize, halfSize);

        float tNear = -std::numeric_limits<float>::max();
        float tFar = std::numeric_limits<float>::max();

        if (fabs(ray.dir.x) < 1e-6) {
            //Ray is parallel to the plane in this axis
            if (ray.start.x < minBounds.x || ray.start.x > maxBounds.x) {
                // Ray is outside the cube, no intersection
                return false;
            }
        } 
        else {
            float t1 = (minBounds.x - ray.start.x) / ray.dir.x;
            float t2 = (maxBounds.x - ray.start.x) / ray.dir.x;

            if (t1 > t2) std::swap(t1, t2);
            if (t1 > tNear) tNear = t1;
            if (t2 < tFar) tFar = t2;

            if (tNear > tFar) return false;
            if (tFar < 0) return false;
        }

        if (fabs(ray.dir.y) < 1e-6) {
            // Ray is parallel to the plane in this axis
            if (ray.start.y < minBounds.y || ray.start.y > maxBounds.y) {
                // Ray is outside the cube, no intersection
                return false;
            }
        } 
        else {
            float t1 = (minBounds.y - ray.start.y) / ray.dir.y;
            float t2 = (maxBounds.y - ray.start.y) / ray.dir.y;

            if (t1 > t2) std::swap(t1, t2);
            if (t1 > tNear) tNear = t1;
            if (t2 < tFar) tFar = t2;

            if (tNear > tFar) return false;
            if (tFar < 0) return false;
        }

        if (fabs(ray.dir.z) < 1e-6) {
            // Ray is parallel to the plane in this axis
            if (ray.start.z < minBounds.z || ray.start.z > maxBounds.z) {
                // Ray is outside the cube, no intersection
                return false;
            }
        } 
        else {
            float t1 = (minBounds.z - ray.start.z) / ray.dir.z;
            float t2 = (maxBounds.z - ray.start.z) / ray.dir.z;

            if (t1 > t2) std::swap(t1, t2);
            if (t1 > tNear) tNear = t1;
            if (t2 < tFar) tFar = t2;

            if (tNear > tFar) return false;
            if (tFar < 0) return false;
        }
        if (tNear > 0) {
            return true;
        } else {
            return false;
        }

    }
    Point pointOfIntersection(Ray ray)
    {
        Point center = bottomLeft + Point(side / 2, side / 2, side / 2);
        double halfSize = side / 2;
        Point minBounds = center - Point(halfSize, halfSize, halfSize);
        Point maxBounds = center + Point(halfSize, halfSize, halfSize);

        float tNear = -std::numeric_limits<float>::max();
        float tFar = std::numeric_limits<float>::max();

        if (fabs(ray.dir.x) < 1e-6) {
            // Ray is parallel to the plane in this axis
            // if (ray.start.x < minBounds.x || ray.start.x > maxBounds.x) {
            //     // Ray is outside the cube, no intersection
            //     return false;
            // }
        } 
        else {
            float t1 = (minBounds.x - ray.start.x) / ray.dir.x;
            float t2 = (maxBounds.x - ray.start.x) / ray.dir.x;

            if (t1 > t2) std::swap(t1, t2);
            if (t1 > tNear) tNear = t1;
            if (t2 < tFar) tFar = t2;

            // if (tNear > tFar) return false;
            // if (tFar < 0) return false;
        }

        if (fabs(ray.dir.y) < 1e-6) {
            // Ray is parallel to the plane in this axis
            // if (ray.start.y < minBounds.y || ray.start.y > maxBounds.y) {
            //     // Ray is outside the cube, no intersection
            //     return false;
            // }
        } 
        else {
            float t1 = (minBounds.y - ray.start.y) / ray.dir.y;
            float t2 = (maxBounds.y - ray.start.y) / ray.dir.y;

            if (t1 > t2) std::swap(t1, t2);
            if (t1 > tNear) tNear = t1;
            if (t2 < tFar) tFar = t2;

            // if (tNear > tFar) return false;
            // if (tFar < 0) return false;
        }

        if (fabs(ray.dir.z) < 1e-6) {
            // Ray is parallel to the plane in this axis
            // if (ray.start.z < minBounds.z || ray.start.z > maxBounds.z) {
            //     // Ray is outside the cube, no intersection
            //     return false;
            // }
        } 
        else {
            float t1 = (minBounds.z - ray.start.z) / ray.dir.z;
            float t2 = (maxBounds.z - ray.start.z) / ray.dir.z;

            if (t1 > t2) std::swap(t1, t2);
            if (t1 > tNear) tNear = t1;
            if (t2 < tFar) tFar = t2;

            // if (tNear > tFar) return false;
            // if (tFar < 0) return false;
        }

        Point intersectionPoint = ray.start + ray.dir * tNear;
        return intersectionPoint;
    }
    
};
vector<Cube> cubes;

class NormalLight
{
public:
    Point position;
    double falloff;
    NormalLight() {};
    NormalLight(Point position, double falloff) : position(position), falloff(falloff) {}
    Point getDirection(Point point)
    {
        return normalize(position - point);
    }
    double getDistance(Point point)
    {
        return distance(position, point);
    }

    bool isPointInShadow(Point point)
    {
        Point dir = normalize(point - position);
        Ray ray(position, dir);
        for(auto sphere : spheres)
        {
            if(sphere.doesIntersect(ray))
            {
                Point intersectionPoint = sphere.pointOfIntersection(ray);
                assert(dcmp(sphere.radius - distance(intersectionPoint, sphere.center)) == 0);
                if(dcmp(distance(intersectionPoint, position) - distance(position, point)) == -1)
                {
                    dbg1(point);
                    dbg1(intersectionPoint);
                    dbg1("sphere in the way");
                    return true;
                }
            }
        }
        for(auto cube : cubes)
        {
            if(cube.doesIntersect(ray))
            {
                Point intersectionPoint = cube.pointOfIntersection(ray);
                if(distance(intersectionPoint, position) + EPS < distance(position, point))
                {
                    return true;
                }
            }
        }
        for(auto pyramid: pyramids)
        {
            if(pyramid.doesIntersect(ray))
            {
                Point intersectionPoint = pyramid.pointOfIntersection(ray);
                if(distance(intersectionPoint, position) + EPS < distance(position, point))
                {
                    return true;
                }
            }
        }
        return false;
    }

    Point reflect(const Point incidentDir, const Point normal) {
        return incidentDir - normal * dot(incidentDir, normal) * 2;
    }

    Point colorAtPoint(Ray ray, Point point, Point normal, Point originalColor, double diffuse, double specular, double shininess)
    {
        Point lightDir = getDirection(point);
        Point reflectionDir = normalize(reflect(ray.dir, normal));
        Point color = Point(0, 0, 0);

        if(isPointInShadow(point)){
            dbg1(ray.object);
            dbg1(spheres.front().pointOfIntersection(ray)); 
            return color;
        } 
        double distance = getDistance(point);
        double scallingFactor = pow(2.71828, -distance * distance * falloff);
        double lambert = max(dot(lightDir, normal) * scallingFactor, 0.0);
        double phong = max(pow(dot(reflectionDir, lightDir), shininess) * scallingFactor, 0.0);

        color = color + originalColor * (lambert * diffuse + phong * specular);

        return color;
    
    
    }
};
vector<NormalLight> lights;

class SpotLight
{
public:
    Point position;
    double falloff;
    Point look;
    double cutOff;
    Point forwadDir;

    SpotLight(){}
    SpotLight(Point position, double fallOff, Point look, double cutOff) : position(position), falloff(falloff), look(look), cutOff(cutOff) 
    {
        forwadDir = normalize(look - position);
    }

    Point getDirection(Point point)
    {
        return normalize(position - point);
    }
    double getDistance(Point point)
    {
        return distance(position, point);
    }

    bool isPointInShadow(Point point)
    {
        Point dir = normalize(point - position);
        Ray ray(position, dir);

        if(acos(dot(dir, forwadDir)) > degreeToRadian(cutOff)) return true;

        for(auto sphere : spheres)
        {
            if(sphere.doesIntersect(ray))
            {
                Point intersectionPoint = sphere.pointOfIntersection(ray);
                assert(dcmp(sphere.radius - distance(intersectionPoint, sphere.center)) == 0);
                if(dcmp(distance(intersectionPoint, position) - distance(position, point)) == -1)
                {
                    dbg1(point);
                    dbg1(intersectionPoint);
                    dbg1("sphere in the way");
                    return true;
                }
            }
        }
        for(auto cube : cubes)
        {
            if(cube.doesIntersect(ray))
            {
                Point intersectionPoint = cube.pointOfIntersection(ray);
                if(distance(intersectionPoint, position) + EPS < distance(position, point))
                {
                    return true;
                }
            }
        }
        for(auto pyramid: pyramids)
        {
            if(pyramid.doesIntersect(ray))
            {
                Point intersectionPoint = pyramid.pointOfIntersection(ray);
                if(distance(intersectionPoint, position) + EPS < distance(position, point))
                {
                    return true;
                }
            }
        }
        return false;
    }

    Point reflect(const Point incidentDir, const Point normal) {
        return incidentDir - normal * dot(incidentDir, normal) * 2;
    }

    Point colorAtPoint(Ray ray, Point point, Point normal, Point originalColor, double diffuse, double specular, double shininess)
    {
        Point lightDir = getDirection(point);
        Point reflectionDir = normalize(reflect(ray.dir, normal));
        Point color = Point(0, 0, 0);

        if(isPointInShadow(point)){
            dbg1(ray.object);
            dbg1(spheres.front().pointOfIntersection(ray)); 
            return color;
        } 
        double distance = getDistance(point);
        double scallingFactor = pow(2.71828, -distance * distance * falloff);
        double lambert = max(dot(lightDir, normal) * scallingFactor, 0.0);
        double phong = max(pow(dot(reflectionDir, lightDir), shininess) * scallingFactor, 0.0);

        color = color + originalColor * (lambert * diffuse + phong * specular);

        return color;
    
    
    }
};
vector<SpotLight> spotLights;


Point rayTrace(Ray ray, int level)
{
    Point color = Point(0, 0, 0);
    Point intersectionPoint;
    Point normal = Point(0, 0, 1);
    
    double t = INT_MAX;

    string object = "none";
    const double ADV = .001;
    if(checkerBoard.doesIntersect(ray))
    {
        intersectionPoint = checkerBoard.intersection(ray);
        ray.object = "checkerboard";
        t = distance(intersectionPoint, ray.start);

        Point reflectDir = normalize(ray.dir - normal * dot(ray.dir, normal) * 2);

        Ray reflectedRay(intersectionPoint + reflectDir * ADV, reflectDir);

        Point originalColor = checkerBoard.getColor(intersectionPoint);
        color = originalColor * checkerBoard.ambient;
        for(NormalLight light: lights)
        {
            color = color + light.colorAtPoint(ray, intersectionPoint, normal, originalColor, checkerBoard.diffuse, 0, 0);
        }
        for(SpotLight spotLight: spotLights)
        {
            color = color + spotLight.colorAtPoint(ray, intersectionPoint, normal, originalColor, checkerBoard.diffuse, 0, 0);
        }
        if(level > 1)
        {
            color = color + rayTrace(reflectedRay, level - 1) * checkerBoard.reflection;
        }
        object = "checkerboard";
    }

    for(auto sphere : spheres)
    {
        if(sphere.doesIntersect(ray))
        {
            Point point = sphere.pointOfIntersection(ray);
            ray.object = "sphere";
            if(distance(point, ray.start) < t)
            {
                t = distance(point, ray.start);
                intersectionPoint = point;
                normal = sphere.normalAtPoint(point);

                Point reflectDir = normalize(ray.dir - normal * dot(ray.dir, normal) * 2);

                Ray reflectedRay(intersectionPoint + reflectDir * ADV, reflectDir);

                Point originalColor = Point(sphere.r, sphere.g, sphere.b);
                color = originalColor * sphere.ambient;
                //dbg1(color);
                for(NormalLight light: lights)
                {
                    assert(dcmp(sphere.radius - distance(intersectionPoint, sphere.center)) == 0);
                    color = color + light.colorAtPoint(ray, intersectionPoint, normal, originalColor, sphere.diffuse, sphere.specular, sphere.shininess);
                    dbg1(color);
                }
                for(SpotLight spotLight: spotLights)
                {
                    assert(dcmp(sphere.radius - distance(intersectionPoint, sphere.center)) == 0);
                    color = color + spotLight.colorAtPoint(ray, intersectionPoint, normal, originalColor, sphere.diffuse, sphere.specular, sphere.shininess);
                    dbg1(color);
                }
                if(level > 1)
                {
                    color = color + rayTrace(reflectedRay, level - 1) * sphere.reflection;
                }
                object = "sphere";
            }
        }
    }
    for(auto cube : cubes)
    {
        if(cube.doesIntersect(ray))
        {
            Point point = cube.pointOfIntersection(ray);
            if(distance(point, ray.start) < t)
            {
                t = distance(point, ray.start);
                intersectionPoint = point;
                normal = cube.normalAtPoint(point);

                Point reflectDir = normalize(ray.dir - normal * dot(ray.dir, normal) * 2);

                Ray reflectedRay(intersectionPoint + reflectDir * ADV, reflectDir);

                Point originalColor = Point(cube.r, cube.g, cube.b);
                color = originalColor * cube.ambient;
                for(NormalLight light: lights)
                {
                    color = color + light.colorAtPoint(ray, intersectionPoint, normal, originalColor, checkerBoard.diffuse, cube.specular, cube.shininess);
                }
                for(SpotLight SpotLight: spotLights)
                {
                    color = color + SpotLight.colorAtPoint(ray, intersectionPoint, normal, originalColor, checkerBoard.diffuse, cube.specular, cube.shininess);
                }
                if(level > 1)
                {
                    color = color + rayTrace(reflectedRay, level - 1) * cube.reflection;
                }
                object = "cube";
            }
        }
    }

    for(Pyramid pyramid: pyramids)
    {
        if(pyramid.doesIntersect(ray))
        {
            Point point = pyramid.pointOfIntersection(ray);
            if(distance(point, ray.start) < t)
            {
                t = distance(point, ray.start);
                intersectionPoint = point;
                normal = pyramid.normalAtPonit(point);

                Point reflectDir = normalize(ray.dir - normal * dot(ray.dir, normal) * 2);

                Ray reflectedRay(intersectionPoint + reflectDir * ADV, reflectDir);

                Point originalColor = Point(pyramid.r, pyramid.g, pyramid.b);
                color = originalColor * pyramid.ambient;
                for(NormalLight light: lights)
                {
                    color = color + light.colorAtPoint(ray, intersectionPoint, normal, originalColor, pyramid.diffuse, pyramid.specular, pyramid.shininess);
                }
                for(SpotLight SpotLight: spotLights)
                {
                    color = color + SpotLight.colorAtPoint(ray, intersectionPoint, normal, originalColor, pyramid.diffuse, pyramid.specular, pyramid.shininess);
                }
                if(level > 1)
                {
                    color = color + rayTrace(reflectedRay, level - 1) * pyramid.reflection;
                }
                object = "pyramid";
            }
        }
    }

    //cout << object << endl;
    if(t > far) return Point(0, 0, 0);
    
    return color;
}


void generateRays()
{
    Point cameraSideDir = normalize(cross(cameraRef - cameraPos, cameraUp));
    Point cameraForwardDir = normalize(cameraRef - cameraPos);
    Point cameraUpDir = normalize(cameraUp);


    Point topLeft = cameraPos + cameraForwardDir * near + cameraUpDir * near * tan(degreeToRadian(fovY / 2)) - cameraSideDir * near * tan(degreeToRadian(fovY * aspectRatio / 2));
    Point bottomLeft = cameraPos + cameraForwardDir * near - cameraUpDir * near * tan(degreeToRadian(fovY / 2)) - cameraSideDir * near * tan(degreeToRadian(fovY * aspectRatio / 2));
    Point topRight = cameraPos + cameraForwardDir * near + cameraUpDir * near * tan(degreeToRadian(fovY / 2)) + cameraSideDir * near * tan(degreeToRadian(fovY * aspectRatio / 2));



    Point rightDir = (topRight - topLeft) / resolution;
    Point bottomDir = (bottomLeft - topLeft) / resolution;
    
    cout << topLeft << endl;
    cout << rightDir << endl;
    cout << bottomDir << endl;

    cout << distance(topLeft, topRight) << endl;
    cout << distance(topLeft, bottomLeft) << endl;

    bitmap_image image(resolution, resolution);
    for (int i = 0; i < resolution; i++)
    {
        for (int j = 0; j < resolution; j++)
        {
            image.set_pixel(i, j, 0, 0, 0);
        }
    }

    for(int i = 0; i < resolution; i++)
    {
        for(int j = 0; j < resolution; j++)
        {
            Point rayStart = topLeft + rightDir * j + bottomDir * i + rightDir / 2 + bottomDir / 2;
            Point rayDir = normalize(rayStart - cameraPos);

            Ray ray(rayStart, rayDir);

            Point color = rayTrace(ray, recursionLevel);
            image.set_pixel(j, i, color.x * 255, color.y * 255, color.z * 255);

        }
    }

    image.save_image("ambient.bmp");
}

void display()
{

    // clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0); // color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    // load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    // initialize the matrix
    glLoadIdentity();

    // now give three info
    // 1. where is the camera (viewer)?
    // 2. where is the camera looking?
    // 3. Which direction is the camera's UP direction?

    gluLookAt(cameraPos.x, cameraPos.y, cameraPos.z, cameraRef.x, cameraRef.y, cameraRef.z, cameraUp.x, cameraUp.y, cameraUp.z);

    // again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);

    /****************************
    / Add your objects from here
    ****************************/
    // add objects
    drawGrid();
    checkerBoard.drawCheckerboard(200);
    for(auto sphere : spheres)
    {
        sphere.drawSphereSegments();
    }
    for(auto cube : cubes)
    {
        cube.DrawCube();
    }
    for(auto pyramid : pyramids)
    {
        pyramid.drawPyramid();
    }
    for(NormalLight light: lights)
    {
        glColor3f(.5, .5, .5);
        glPushMatrix();
        glTranslatef(light.position.x, light.position.y, light.position.z);
        glutSolidSphere(10, 100, 100);
        glPopMatrix();
    }
    for(SpotLight spotLight: spotLights)
    {
        Point curDir = Point(0, 0, 1);
        double angle = -acos(dot(normalize(spotLight.position - spotLight.look), curDir)) * 180 / PI;
        Point axis = normalize(cross(spotLight.position - spotLight.look, curDir));

        glColor3d(1, 1, 0);
        glPushMatrix();
        glTranslated(spotLight.position.x, spotLight.position.y, spotLight.position.z); 
        glRotated(angle, axis.x, axis.y, axis.z);
        glutSolidCone(20, 20 / tan(degreeToRadian(spotLight.cutOff)), 100, 100);
        glPopMatrix();
    }
    // ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}

void animate()
{
    // codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init()
{
    // codes for initialization
    cameraPos = Point(-100, -150, 10);
    cameraRef = Point(0, 1, 0);
    cameraUp = Point(0, 0, 1);
    angle = 0;

    // clear the screen
    glClearColor(0, 0, 0, 0);

    /************************
    / set-up projection here
    ************************/
    // load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    // initialize the matrix
    glLoadIdentity();

    // give PERSPECTIVE parameters
    gluPerspective(fovY, aspectRatio, near, far);
    // field of view in the Y (vertically)
    // aspect ratio that determines the field of view in the X direction (horizontally)
    // near distance
    // far distance

}

void takeInput()
{
    ifstream fin("description.txt");
    fin >> near >> far;
    fin >> fovY >> aspectRatio;
    fin >> recursionLevel;
    fin >> resolution;

    fin >> checkerBoard.cellSize;
    fin >> checkerBoard.ambient >> checkerBoard.diffuse >> checkerBoard.reflection;

    fin >> objectCount;

    while(objectCount--)
    {
        string object;
        fin >> object;

        if(object == "sphere")
        {
            Point center;
            double radius;
            double r, g, b;
            double ambient, diffuse, specular, reflection, shininess;
            fin >> center.x >> center.y >> center.z >> radius;
            fin >> r >> g >> b;
            dbg1(r, g, b);
            fin >> ambient >> diffuse >> specular >> reflection >> shininess;
            spheres.push_back(Spherre(center, radius, r, g, b, ambient, diffuse, specular, reflection, shininess));
        }
        else if(object == "cube")
        {
            Point bottomLeft;
            double side;
            double r, g, b;
            double ambient, diffuse, specular, reflection, shininess;
            fin >> bottomLeft.x >> bottomLeft.y >> bottomLeft.z >> side;
            fin >> r >> g >> b;
            fin >> ambient >> diffuse >> specular >> reflection >> shininess;
            cubes.push_back(Cube(bottomLeft, side, r, g, b, ambient, diffuse, specular, reflection, shininess));
        }
        else if(object == "pyramid")
        {
            Point bottomLeft;
            double side;
            double height;
            double r, g, b;
            double ambient, diffuse, specular, reflection, shininess;
            fin >> bottomLeft.x >> bottomLeft.y >> bottomLeft.z >> side >> height;
            fin >> r >> g >> b;
            fin >> ambient >> diffuse >> specular >> reflection >> shininess;
            //pyramids.push_back(Pyramid(bottomLeft, side, height, r, g, b, ambient, diffuse, specular, reflection, shininess));
        }
    }

    int normalLightCount;
    fin >> normalLightCount;

    while(normalLightCount--)
    {
        double normalLightX, normalLightY, normalLightZ;
        double falloff;
        fin >> normalLightX >> normalLightY >> normalLightZ;
        fin >> falloff;

        

        lights.push_back(NormalLight(Point(normalLightX, normalLightY, normalLightZ), falloff));

    }

    int spotLightCount;
    fin >> spotLightCount;

    while(spotLightCount--)
    {
        double spotLightX, spotLightY, spotLightZ;
        double falloff;
        Point look;
        double cutOff;
        fin >> spotLightX >> spotLightY >> spotLightZ;
        fin >> falloff;
        fin >> look.x >> look.y >> look.z;
        fin >> cutOff;
        

        spotLights.push_back(SpotLight(Point(spotLightX, spotLightY, spotLightZ), falloff, look, cutOff));

    }
}

int main(int argc, char **argv)
{
    takeInput();
    glutInit(&argc, argv);
    glutInitWindowSize(resolution, resolution);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // Depth, Double buffer, RGB color

    glutCreateWindow("Magic Cube");

    init();

    glEnable(GL_DEPTH_TEST); // enable Depth Testing

    glutDisplayFunc(display); // display callback function
    glutIdleFunc(animate);    // what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop(); // The main loop of OpenGL

    return 0;
}
