/*
CSCI 420
Assignment 3 Raytracer

Name: Junmeng Xu
*/
#include <stdlib.h>
#include <string.h>
#include "pic.h"
#include <math.h>
#include <limits.h>

// For Linux
//#include <GL/gl.h>
//#include <GL/glu.h>
//#include <GL/glut.h>

// For Mac
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//For debugging
//#define WIDTH 200
//#define HEIGHT 150

//the field of view of the camera
#define fov 60.0

//the math constant
#define PI 3.14159265

//For generate sphere light from point light
#define POINTS_FROM_SPHERE 8
#define DISTANCE_FROM_CENTER 0.1

//For math calculation
double RADIAN = (double) PI/180;
double EPSILON = pow(10, -14);
double MAX_Z  = -pow(10, 14);

//Times of reflections
double numReflection = 0;

//Enable rendering soft shadows
bool g_SoftShadows = false;

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

/*
 * Represent vector mostly or point sometimes
 */
struct Vector {
   double x;
   double y;
   double z;
    
    /*
     * Define some basic functions for vector calculation
     */
    
    /* cross product for vectors */
    Vector vectorCross(Vector v2)
    {
        Vector v;
        v.x = y*v2.z - z*v2.y;
        v.y = z*v2.x - x*v2.z;
        v.z = x*v2.y - y*v2.x;
        return v;
    }

    /* dot product for vectors */
    double vectorDot(Vector v2)
    {
        double ans = x*v2.x + y*v2.y + z*v2.z;
        return ans;
    }
    
    /* calculate the normal vector of input vector */
    Vector vectorNormal()
    {
        Vector v;
        double len = sqrt(x*x + y*y + z*z);
        v.x = x/len;
        v.y = y/len;
        v.z = z/len;
        return v;
    }

    /* calculate the length of one vector */
    double vectorLength()
    {
        double len = sqrt(x*x + y*y + z*z);
        return len;
    }
    
    /* vector addition */
    Vector vectorAdd(Vector vector)
    {
        Vector v;
        v.x = x + vector.x;
        v.y = y + vector.y;
        v.z = z + vector.z;
        return v;
    }
    
    /* vector subtraction */
    Vector vectorSub(Vector vector)
    {
        Vector v;
        v.x = x - vector.x;
        v.y = y - vector.y;
        v.z = z - vector.z;
        return v;
    }
    
    /* vector scaler */
    Vector vectorScaler(double s)
    {
        Vector v;
        v.x = s * x;
        v.y = s * y;
        v.z = s * z;
        return v;
    }
    
    /* vector reverse */
    Vector vectorReverse()
    {
        Vector v;
        v.x = -x;
        v.y = -y;
        v.z = -z;
        return v;
    }
};

/*
 * Represent color with r g b
 */
struct Color {
    double r;
    double g;
    double b;
    
    /*
     Define some basic functions for color and light calculation
     */
    
    /* color add a light */
    Color addLight(double light[3])
    {
        Color color;
        color.r = r + light[0];
        color.g = g + light[1];
        color.b = b + light[2];
        return color;
    }
    
    /* color add another color */
    Color addColor(Color c)
    {
        Color color;
        color.r = r + c.r;
        color.g = g + c.g;
        color.b = b + c.b;
        return color;
    }
    
    /* Clamp color's number when it is more than 1 or less than 0 */
    void clamp()
    {
        if(r > 1.0f){
            r = 1.0f;
        }else if(r < 0.0f){
            r = 0.0f;
        }
        
        if(g > 1.0f){
            g = 1.0f;
        }else if(g < 0.0f){
            r = 0.0f;
        }
        
        if(b > 1.0f){
            b = 1.0f;
        }else if(b < 0.0f){
            b = 0.0f;
        }
    }
};

Color DARK = {0.0f, 0.0f, 0.0f};
Color WHITE = {1.0f, 1.0f, 1.0f};

typedef struct _Sphere_Light
{
  double position[3];
  double color[3];
  Vector dir;
} SphereLight;

/*
 * Represent Ray tracing
 * contains origin and direction
 */
struct rayTracer
{
    Vector origin;
    Vector direction;
    
    /* Check the intersection between the ray and sphere */
    bool checkSphereIntersection(Sphere sphere, Vector& point)
    {
        double r = sphere.radius;
        Vector vCtoO;
        vCtoO.x = origin.x - sphere.position[0];
        vCtoO.y = origin.y - sphere.position[1];
        vCtoO.z = origin.z - sphere.position[2];
        
        double a = 1;
        double b = 2 * direction.vectorDot(vCtoO);
        double c = vCtoO.vectorDot(vCtoO) - pow(r, 2);
        
        double numeratorRight = b*b - 4*a*c;
        double t0, t1;
        t0 = MAX_Z;
        t1 = MAX_Z;
        
        // no real solution for the equation
        if(numeratorRight < 0)
        {
            return false;
        }
        
        double absNumeratorRight = fabs(numeratorRight);

        t0 = (-b + sqrt(absNumeratorRight)) * 0.5f;
        t1 = (-b - sqrt(absNumeratorRight)) * 0.5f;
        
        // no intersection
        if(t0 < 0 && t1 < 0) return false;
        
        // find the nearest intersection
        double t = t0;
        if(t1 > 0){
            t = t1;
        }
        point = origin.vectorAdd(direction.vectorScaler(t));
        return true;
    }
    
    /* Check the intersection between the ray and triangle
       Intersect ray with plane containing polygon
       Check if intersection point is inside polygon
     */
    bool checkTriangleIntersection(Triangle triangle, Vector& point)
    {
        Vector A = {triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]};
        Vector B = {triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]};
        Vector C = {triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]};
        
        Vector N;
        N = B.vectorSub(A).vectorCross(C.vectorSub(A)).vectorNormal();
        double NdotD = N.vectorDot(direction);
        
        // no intersection (ray parallel to plane)
        if(fabs(NdotD) < EPSILON) return false;
        
        double d = - N.vectorDot(A);
        double t = - (N.vectorDot(origin) + d) / NdotD;
        
        // the intersection is behind ray origin
        if(t <= 0) return false;
        
        point = origin.vectorAdd( direction.vectorScaler(t) );
        
        double alpha = B.vectorSub(A).vectorCross(point.vectorSub(A)).vectorDot(N);
        double beta = C.vectorSub(B).vectorCross(point.vectorSub(B)).vectorDot(N);
        double gamma = A.vectorSub(C).vectorCross(point.vectorSub(C)).vectorDot(N);
        
        // intersection point is outside of triangle
        if(alpha < 0 || beta < 0 || gamma < 0) return false;
        
        return true;
    }
};

/*
 Uniformly send out rays from the camera location center
 of projection (COP)
 */
rayTracer fireRays(double pixelX, double pixelY)
{
    double aspectRatio = (double) WIDTH/HEIGHT;
    double tanFovHalf = tan((fov * 0.5f) * RADIAN);
    
    // the camera is in 0,0,0
    // screen coordinates are [-1,1]
    // so we need to map x,y into this area
    rayTracer primaryRay;
    pixelX = (pixelX + 0.1) / WIDTH;
    pixelY = (pixelY + 0.1) / HEIGHT;
    double pixelXScreen = 2.0f * pixelX - 1;
    double pixelYScreen = 2.0f * pixelY - 1;
    
    // map screen coordinates into image plan
    double pixelXCamera = pixelXScreen * tanFovHalf * aspectRatio;
    double pixelYCamera = pixelYScreen * tanFovHalf;
    
    // get ray vector
    Vector rayVector;
    rayVector.x = pixelXCamera;
    rayVector.y = pixelYCamera;
    rayVector.z =  -1.0;
    
    // set ray origin as [0,0,0]
    Vector origin;
    origin.x = 0;
    origin.y = 0;
    origin.z = 0;
    
    primaryRay.origin = origin;
    primaryRay.direction = rayVector.vectorNormal();
    
    return primaryRay;
}

/*
 Clamp color
 if it is greater than 1.0, clamp it to 1.0
 if it is less than 0.0, clamp it to 0.0
 Color values range from 0-1
 */
double clamp(double in)
{
    if(in > 1){
        return 1.0;
    }else if(in < 0){
        return 0.0;
    }
    return in;
}

/*
 Phong shading for sphere
 I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ α)
 */
Color spherePhongShading(Sphere sphere, Light light, Vector intersection)
{
    Vector light_position;
    light_position.x = light.position[0];
    light_position.y = light.position[1];
    light_position.z = light.position[2];
    
    Vector sphere_position;
    sphere_position.x = sphere.position[0];
    sphere_position.y = sphere.position[1];
    sphere_position.z = sphere.position[2];
    
    Vector L = light_position.vectorSub(intersection).vectorNormal();
    Vector N = intersection.vectorSub(sphere_position).vectorNormal();
    double LdotN = L.vectorDot(N);
    // L dot N, or R dot V are negative, clamp them to zero
    LdotN = clamp(LdotN);
    
    Vector R = (N.vectorScaler( 2 * (LdotN) ).vectorSub(L) ).vectorNormal();
    Vector V = intersection.vectorReverse().vectorNormal();
    double RdotV = R.vectorDot(V);
    // L dot N, or R dot V are negative, clamp them to zero
    RdotV = clamp(RdotV);
    
    // for each color channel separately
    
    // the diffuse reflection coefficient kd
    Color kd;
    kd.r = sphere.color_diffuse[0];
    kd.g = sphere.color_diffuse[1];
    kd.b = sphere.color_diffuse[2];
    
    // the Specular reflection coefficient ks
    Color ks;
    ks.r = sphere.color_specular[0];
    ks.g = sphere.color_specular[1];
    ks.b = sphere.color_specular[2];
    
    double specular = pow(RdotV, sphere.shininess);
    
    Color l;
    l.r = light.color[0] * (kd.r * LdotN + ks.r * specular);
    l.g = light.color[1] * (kd.g * LdotN + ks.g * specular);
    l.b = light.color[2] * (kd.b * LdotN + ks.b * specular);
    
    return l;
}

/*
 Phong shading for triangle
 I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ α)
 interpolate the vertex normals given inscene file to the specific location of
 the ray-triangle intersection (using barycentric coordinates)
 */
Color trianglePhongShading(Triangle triangle, Light light, Vector intersection)
{
    Vector A = {triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]};
    Vector B = {triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]};
    Vector C = {triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]};
    
    // using barycentric coordinates
    double triangleArea = B.vectorSub(A).vectorCross( C.vectorSub(A) ).vectorLength();
    double alpha = B.vectorSub(A).vectorCross( intersection.vectorSub(A) ).vectorLength() / triangleArea;
    double beta = C.vectorSub(B).vectorCross( intersection.vectorSub(B) ).vectorLength() / triangleArea;
    double gamma = 1.0-alpha-beta;
    
    // interpolate the vertex normals
    Vector N;
    N.x = alpha * triangle.v[2].normal[0] + beta * triangle.v[0].normal[0] + gamma * triangle.v[1].normal[0];
    N.y = alpha * triangle.v[2].normal[1] + beta * triangle.v[0].normal[1] + gamma * triangle.v[1].normal[1];
    N.z = alpha * triangle.v[2].normal[2] + beta * triangle.v[0].normal[2] + gamma * triangle.v[1].normal[2];
    N = N.vectorNormal();
    
    // unit vector to light
    Vector light_position = {light.position[0], light.position[1], light.position[2]};
    Vector L = light_position.vectorSub(intersection).vectorNormal();
    double LdotN = L.vectorDot(N);
    LdotN = clamp(LdotN);
    
    // L reflected about N
    Vector R = N.vectorScaler(2 * LdotN).vectorSub(L).vectorNormal();
    Vector V = intersection.vectorReverse().vectorNormal();
    double RdotV = R.vectorDot(V);
    RdotV = clamp(RdotV);
    
    // interpolate the diffuse reflection coefficient kd (material)
    Color kd;
    kd.r = alpha * triangle.v[2].color_diffuse[0] + beta * triangle.v[0].color_diffuse[0] + gamma * triangle.v[1].color_diffuse[0];
    kd.g = alpha * triangle.v[2].color_diffuse[1] + beta * triangle.v[0].color_diffuse[1] + gamma * triangle.v[1].color_diffuse[1];
    kd.b = alpha * triangle.v[2].color_diffuse[2] + beta * triangle.v[0].color_diffuse[2] + gamma * triangle.v[1].color_diffuse[2];
    
    // interpolate the Specular reflection coefficient ks (material)
    Color ks;
    ks.r = alpha * triangle.v[2].color_specular[0] + beta * triangle.v[0].color_specular[0] + gamma * triangle.v[1].color_specular[0];
    ks.g = alpha * triangle.v[2].color_specular[1] + beta * triangle.v[0].color_specular[1] + gamma * triangle.v[1].color_specular[1];
    ks.b = alpha * triangle.v[2].color_specular[2] + beta * triangle.v[0].color_specular[2] + gamma * triangle.v[1].color_specular[2];
    
    // interpolate the shininess coefficient
    double shininess = alpha * triangle.v[2].shininess + beta * triangle.v[0].shininess + gamma * triangle.v[1].shininess;
    double specular = pow(RdotV, shininess);
    
    // get the final color
    Color l;
    l.r = light.color[0] * (kd.r * LdotN + ks.r * specular);
    l.g = light.color[1] * (kd.g * LdotN + ks.g * specular);
    l.b = light.color[2] * (kd.b * LdotN + ks.b * specular);
    
    return l;
}

/*
 For each ray, check all intersections among it with all spheres
 */
Color checkAllSphereIntersection(rayTracer ray, Color color, Vector& intersection, double& nearestZ, int& sphereIntersect)
{
    Color pixel_color = color;
    
    // check every sphere
    for(int i=0; i<num_spheres; i++){
        // initial nearest z coordinate
        Vector intersect = {0.0, 0.0, MAX_Z};
        
        //  if the intersection's z with sphere is more near than previous one
        if(ray.checkSphereIntersection(spheres[i], intersect) && intersect.z > nearestZ){
            pixel_color = DARK;
            sphereIntersect = i;
            
            // check every light
            for(int j=0; j<num_lights; j++){
                Vector light;
                light.x = lights[j].position[0];
                light.y = lights[j].position[1];
                light.z = lights[j].position[2];
                
                // set the intersecton point as the origin
                // check whether the point is in the shadow
                // follow the direction to check intersections with other objects
                Vector shadowDirection = light.vectorSub(intersect).vectorNormal();
                rayTracer shadow;
                shadow.origin = intersect;
                shadow.direction = shadowDirection;
                
                bool inShadow = false;
                
                // check whether the ray intersects with sphere
                for(int k=0; k<num_spheres; k++){
                    Vector hitSphere;
                    if(shadow.checkSphereIntersection(spheres[k], hitSphere) && (k != i))
                    {
                        Vector shadowFromOriginToHitSphere = hitSphere.vectorSub(intersect);
                        Vector shadowFromOriginToLight = light.vectorSub(intersect);
                        if(shadowFromOriginToLight.vectorLength() > shadowFromOriginToHitSphere.vectorLength())
                        {
                            inShadow = true;
                            break;
                        }
                    }
                }
                
                // check whether the ray intersects with triangle
                for(int k=0; k<num_triangles; k++){
                    Vector hitTriangles;
                    if(shadow.checkTriangleIntersection(triangles[k], hitTriangles)){
                        Vector shadowFromOriginToTriangle = hitTriangles.vectorSub(intersect);
                        Vector shadowFromOriginToLight = light.vectorSub(intersect);
                        if(shadowFromOriginToLight.vectorLength() > shadowFromOriginToTriangle.vectorLength()){
                            inShadow = true;
                            break;
                        }
                    }
                }
                
                // if ray can get to light without intersection
                if(!inShadow){
                    Color temp = pixel_color;
                    pixel_color = temp.addColor(spherePhongShading(spheres[i], lights[j], intersect));
                }
            }
            
            // refresh nearest z and intersection point
            nearestZ = intersect.z;
            intersection = intersect;
        }
    }
    return pixel_color;
}

/*
 For each ray, check all intersections among it with all triangles
 */
Color checkAllTriangleIntersection(rayTracer& ray, const Color& color, Vector& intersection, double& nearestZ, int& triangleIntersect)
{
    Color pixel_color = color;
    
    // check every triangle
    for(int i=0; i<num_triangles; i++){
        // initial nearest z coordinate
        Vector intersect = {0.0, 0.0, MAX_Z};
        
        //  if the intersection's z with sphere is more near than previous one
        if(ray.checkTriangleIntersection(triangles[i], intersect) && intersect.z > nearestZ){
            pixel_color = DARK;
            triangleIntersect = i;
            
            // check every light
            for(int j=0; j<num_lights; j++){
                Vector light;
                light.x = lights[j].position[0];
                light.y = lights[j].position[1];
                light.z = lights[j].position[2];
                
                // set the intersecton point as the origin
                // check whether the point is in the shadow
                // follow the direction to check intersections with other objects
                Vector shadowDirection = light.vectorSub(intersect).vectorNormal();
                rayTracer shadow;
                shadow.origin = intersect;
                shadow.direction = shadowDirection;
                
                bool inShadow = false;
                
                // check whether the ray intersects with sphere
                for(int k=0; k<num_spheres; k++){
                    Vector hitSphere;
                    if(shadow.checkSphereIntersection(spheres[k], hitSphere))
                    {
                        Vector shadowFromOriginToHitSphere = hitSphere.vectorSub(intersect);
                        Vector shadowFromOriginToLight = light.vectorSub(intersect);
                        if(shadowFromOriginToLight.vectorLength() > shadowFromOriginToHitSphere.vectorLength())
                        {
                            inShadow = true;
                            break;
                        }
                    }
                }
                
                // check whether the ray intersects with triangle
                for(int k=0; k<num_triangles; k++){
                    Vector hitTriangles;
                    if(shadow.checkTriangleIntersection(triangles[k], hitTriangles) && k != i){
                        Vector shadowFromOriginToTriangle = hitTriangles.vectorSub(intersect);
                        Vector shadowFromOriginToLight = light.vectorSub(intersect);
                        if(shadowFromOriginToLight.vectorLength() > shadowFromOriginToTriangle.vectorLength()){
                            inShadow = true;
                            break;
                        }
                    }
                }
                
                // if ray can get to light without intersection
                if(!inShadow){
                    Color temp = pixel_color;
                    pixel_color = temp.addColor(trianglePhongShading(triangles[i], lights[j], intersect));
                }
            }
            
            // refresh nearest z and intersection point
            nearestZ = intersect.z;
            intersection = intersect;
        }
    }
    return pixel_color;
}


/*
 *** Extra Credit : Soft Shadows ***
 Generate sphere light for each point light
 In real life, there is no ideal point light
 Also we need sphere light to get soft shadow
 */
Light * transferToSphereLights(Light light)
{
    // Generate [POINTS_FROM_SPHERE*POINTS_FROM_SPHERE] number of point light to mimic
    // sphere light for every original point light
    Light * sphereLights = new Light[POINTS_FROM_SPHERE*POINTS_FROM_SPHERE];
    int index = 0;
    
    // For 3D sphere points generating
    for(int i=0; i<POINTS_FROM_SPHERE; i++){
        double phi = i * (360 / (POINTS_FROM_SPHERE)) * RADIAN;
        double sinPhi = sin(phi);
        double cosPhi = cos(phi);
        for(int j=0; j<POINTS_FROM_SPHERE; j++){
            double theta = j * (360 / (POINTS_FROM_SPHERE)) * RADIAN;
            double sinTheta = sin(theta);
            double cosTheta = cos(theta);

            double x = DISTANCE_FROM_CENTER * sinPhi * cosTheta;
            double y = DISTANCE_FROM_CENTER * sinPhi * sinTheta;
            double z = DISTANCE_FROM_CENTER * cosPhi;

            // determine sphere point's location
            sphereLights[index].position[0] = light.position[0] + x;
            sphereLights[index].position[1] = light.position[1] + y;
            sphereLights[index].position[2] = light.position[2] + z;

            // determine its light color
            // and dilute light strength to prevent bright color
            sphereLights[index].color[0] = light.color[0]/(POINTS_FROM_SPHERE*POINTS_FROM_SPHERE);
            sphereLights[index].color[1] = light.color[1]/(POINTS_FROM_SPHERE*POINTS_FROM_SPHERE);
            sphereLights[index].color[2] = light.color[2]/(POINTS_FROM_SPHERE*POINTS_FROM_SPHERE);

            index++;
        }
    }
    return sphereLights;
}

/*
 *** Extra Credit : Soft Shadows ***
 For each ray, check all intersections among it with all spheres
 here for soft shadows, we need to check sphere light for each point light
 */
Color checkAllSphereIntersectionSoftShadows(rayTracer& ray, const Color& color, Vector& intersection, double& nearestZ, int& sphereIntersect)
{
    Color pixel_color = color;
    
    // check every sphere
    for(int i=0; i<num_spheres; i++){
        Vector intersect = {0.0, 0.0, MAX_Z};
        
        // if the intersection's z with sphere is more near than previous one
        if(ray.checkSphereIntersection(spheres[i], intersect) && intersect.z > nearestZ){
            pixel_color = DARK;
            sphereIntersect = i;
            
            // check every light
            for(int j=0; j<num_lights; j++){
                
                // check each point light from the sphere light generated from original point light
                Light * sphereLights = transferToSphereLights(lights[j]);
                for(int t=0; t<pow(POINTS_FROM_SPHERE,2); t++){
                    // get current light attributes
                    Vector light;
                    light.x = sphereLights[t].position[0];
                    light.y = sphereLights[t].position[1];
                    light.z = sphereLights[t].position[2];
                    
                    Light curLight;
                    curLight.position[0] = sphereLights[t].position[0];
                    curLight.position[1] = sphereLights[t].position[1];
                    curLight.position[2] = sphereLights[t].position[2];
                    curLight.color[0] = sphereLights[t].color[0];
                    curLight.color[1] = sphereLights[t].color[1];
                    curLight.color[2] = sphereLights[t].color[2];
                    
                    // set the intersecton point as the origin
                    // check whether the point is in the shadow
                    // follow the direction to check intersections with other objects
                    Vector shadowDirection = light.vectorSub(intersect).vectorNormal();
                    rayTracer shadow;
                    shadow.origin = intersect;
                    shadow.direction = shadowDirection;
                    
                    bool inShadow = false;
                    
                    // check whether the ray intersects with sphere
                    for(int k=0; k<num_spheres; k++){
                        Vector hitSphere;
                        if(shadow.checkSphereIntersection(spheres[k], hitSphere) && (k != i))
                        {
                            Vector shadowFromOriginToHitSphere = hitSphere.vectorSub(intersect);
                            Vector shadowFromOriginToLight = light.vectorSub(intersect);
                            if(shadowFromOriginToLight.vectorLength() > shadowFromOriginToHitSphere.vectorLength())
                            {
                                inShadow = true;
                                break;
                            }
                        }
                    }
                    
                    // check whether the ray intersects with triangle
                    for(int k=0; k<num_triangles; k++){
                        Vector hitTriangles;
                        if(shadow.checkTriangleIntersection(triangles[k], hitTriangles)){
                            Vector shadowFromOriginToTriangle = hitTriangles.vectorSub(intersect);
                            Vector shadowFromOriginToLight = light.vectorSub(intersect);
                            if(shadowFromOriginToLight.vectorLength() > shadowFromOriginToTriangle.vectorLength()){
                                inShadow = true;
                                break;
                            }
                        }
                    }
                    
                    // if ray can get to light without intersection
                    if(!inShadow){
                        Color temp = pixel_color;
                        pixel_color = temp.addColor(spherePhongShading(spheres[i], curLight, intersect));
                    }
                }
            }
            
            // refresh nearest z and intersection point
            nearestZ = intersect.z;
            intersection = intersect;
        }
    }
    return pixel_color;
}


/*
 *** Extra Credit : Soft Shadows ***
 For each ray, check all intersections among it with all triangles
 here for soft shadows, we need to check sphere light for each point light
 */
Color checkAllTriangleIntersectionSoftShadows(rayTracer& ray, const Color& color, Vector& intersection, double& nearestZ, int& triangleIntersect)
{
    Color pixel_color = color;
    
    // check every triangle
    for(int i=0; i<num_triangles; i++){
        Vector intersect = {0.0, 0.0, MAX_Z};
        
        // if the intersection's z with sphere is more near than previous one
        if(ray.checkTriangleIntersection(triangles[i], intersect) && intersect.z > nearestZ){
            pixel_color = DARK;
            triangleIntersect = i;
            
            // check every light
            for(int j=0; j<num_lights; j++){
                
                // check each point light from the sphere light generated from original point light
                Light * sphereLights = transferToSphereLights(lights[j]);
                for(int t=0; t<pow(POINTS_FROM_SPHERE,2); t++){
                    // get current light attributes
                    Vector light;
                    light.x = sphereLights[t].position[0];
                    light.y = sphereLights[t].position[1];
                    light.z = sphereLights[t].position[2];
                    
                    Light curLight;
                    curLight.position[0] = sphereLights[t].position[0];
                    curLight.position[1] = sphereLights[t].position[1];
                    curLight.position[2] = sphereLights[t].position[2];
                    curLight.color[0] = sphereLights[t].color[0];
                    curLight.color[1] = sphereLights[t].color[1];
                    curLight.color[2] = sphereLights[t].color[2];
                    
                    // set the intersecton point as the origin
                    // check whether the point is in the shadow
                    // follow the direction to check intersections with other objects
                    Vector shadowDirection = light.vectorSub(intersect).vectorNormal();
                    rayTracer shadow;
                    shadow.origin = intersect;
                    shadow.direction = shadowDirection;

                    bool inShadow = false;
                    
                    // check whether the ray intersects with sphere
                    for(int k=0; k<num_spheres; k++){
                        Vector hitSphere;
                        if(shadow.checkSphereIntersection(spheres[k], hitSphere))
                        {
                            Vector shadowFromOriginToHitSphere = hitSphere.vectorSub(intersect);
                            Vector shadowFromOriginToLight = light.vectorSub(intersect);
                            if(shadowFromOriginToLight.vectorLength() > shadowFromOriginToHitSphere.vectorLength())
                            {
                                inShadow = true;
                                break;
                            }
                        }
                    }
                    
                    // check whether the ray intersects with triangle
                    for(int k=0; k<num_triangles; k++){
                        Vector hitTriangles;
                        if(shadow.checkTriangleIntersection(triangles[k], hitTriangles) && k != i){
                            Vector shadowFromOriginToTriangle = hitTriangles.vectorSub(intersect);
                            Vector shadowFromOriginToLight = light.vectorSub(intersect);
                            if(shadowFromOriginToLight.vectorLength() > shadowFromOriginToTriangle.vectorLength()){
                                inShadow = true;
                                break;
                            }
                        }
                    }
                    
                    // if ray can get to light without intersection
                    if(!inShadow){
                        Color temp = pixel_color;
                        pixel_color = temp.addColor(trianglePhongShading(triangles[i], curLight, intersect));
                    }
                }
            }
            
            // refresh nearest z and intersection point
            nearestZ = intersect.z;
            intersection = intersect;
        }
    }
    return pixel_color;
}

/*
 *** Extra Credit : Recursive reflection ***
 call ray tracer recursively
 fire the shadow rays and compute a local Phong color
 The final color should equal (1 - ks) * localPhongColor + ks * colorOfReflectedRay.
 */
Color recursiveReflection(rayTracer& ray, int reflectionLevel)
{
    // if real times of reflection is more than we need, return DARK color
    if(reflectionLevel > numReflection){
        return DARK;
    }
    
    // pre-define variables for later use
    Vector intersection;
    Color fromLight;
    fromLight = WHITE;
    double nearestZ = MAX_Z;
    
    // find the object the ray intersects with refreshing nearest Z and intersection point
    int sphereIndex=-1, triangleIndex=-1;
    if(g_SoftShadows){
        fromLight = checkAllSphereIntersectionSoftShadows(ray, fromLight, intersection, nearestZ, sphereIndex);
        fromLight = checkAllTriangleIntersectionSoftShadows(ray, fromLight, intersection, nearestZ, triangleIndex);
    }else{
        fromLight = checkAllSphereIntersection(ray, fromLight, intersection, nearestZ, sphereIndex);
        fromLight = checkAllTriangleIntersection(ray, fromLight, intersection, nearestZ, triangleIndex);
    }
    
    // find the reflection ray and corresponding
    Color ks;
    rayTracer reflection;
    
    // check whether the reflection ray comes from intersection on sphere
    if(sphereIndex != -1){
        Sphere sphere = spheres[sphereIndex];
        Vector spherePos = {sphere.position[0], sphere.position[1], sphere.position[2]};
        
        // calculate normal and direction
        Vector N = intersection.vectorSub(spherePos).vectorNormal();
        Vector L = {-ray.direction.x, -ray.direction.y, -ray.direction.z};
        double LdotN = L.vectorDot(N);
        LdotN = clamp(LdotN);
        Vector R = N.vectorScaler(2 * LdotN).vectorSub(L).vectorNormal();
        Vector origin = intersection.vectorAdd(R);
        
        // set reflection ray and ks
        reflection.origin = origin;
        reflection.direction = R.vectorNormal();
        ks.r = sphere.color_specular[0];
        ks.g = sphere.color_specular[1];
        ks.b = sphere.color_specular[2];
    }
    
    // check whether the reflection ray comes from intersection on triangle
    // notice if triangleIndex != -1, means the nearest Z for spheres is still more far than
    // its for triangles, so there is a triangle in the more near location
    if(triangleIndex != -1){
        Triangle triangle = triangles[triangleIndex];
        Vector A = {triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]};
        Vector B = {triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]};
        Vector C = {triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]};
        
        // calculate normal and direction
        double triangleArea = B.vectorSub(A).vectorCross( C.vectorSub(A) ).vectorLength();
        double alpha = B.vectorSub(A).vectorCross( intersection.vectorSub(A) ).vectorLength() / triangleArea;
        double beta = C.vectorSub(B).vectorCross( intersection.vectorSub(B) ).vectorLength() / triangleArea;
        double gamma = 1.0-alpha-beta;
        
        Vector N;
        N.x = alpha * triangle.v[2].normal[0] + beta * triangle.v[0].normal[0] + gamma * triangle.v[1].normal[0];
        N.y = alpha * triangle.v[2].normal[1] + beta * triangle.v[0].normal[1] + gamma * triangle.v[1].normal[1];
        N.z = alpha * triangle.v[2].normal[2] + beta * triangle.v[0].normal[2] + gamma * triangle.v[1].normal[2];
        N = N.vectorNormal();
        
        Vector L = {-ray.direction.x, -ray.direction.y, -ray.direction.z};
        double LdotN = L.vectorDot(N);
        LdotN = clamp(LdotN);
        
        Vector R = N.vectorScaler(2 * LdotN).vectorSub(L).vectorNormal();
        Vector origin = intersection.vectorAdd(R.vectorScaler(EPSILON));
        
        // set reflection ray and ks
        reflection.origin = origin;
        reflection.direction = R.vectorNormal();
        
        ks.r = alpha * triangle.v[2].color_specular[0] + beta * triangle.v[0].color_specular[0] + gamma * triangle.v[1].color_specular[0];
        ks.g = alpha * triangle.v[2].color_specular[1] + beta * triangle.v[0].color_specular[1] + gamma * triangle.v[1].color_specular[1];
        ks.b = alpha * triangle.v[2].color_specular[2] + beta * triangle.v[0].color_specular[2] + gamma * triangle.v[1].color_specular[2];
    }
    
    Color fromThisLevel;
    
    // if there is no reflectino, which means only primary ray
    if(numReflection == 0)
    {
        // no intersection with any objects
        if(sphereIndex == -1 && triangleIndex == -1){
            // return default color
            return WHITE;
            
        // there is intersection
        }else{
            // return phong shading color
            return fromLight;
        }
    }
    // there are reflections as well as primary ray
    else
    {
        // no intersection with any objects
        if(sphereIndex == -1 && triangleIndex == -1){
            // for primary ray, return default color
            if(reflectionLevel == 0) return WHITE;
            // for reflection ray, return DARK, that is no color
            return DARK;
        }
        // for reflection
        // the final color should equal (1 - ks) * localPhongColor + ks * colorOfReflectedRay.
        else{
            Color fromReflection = recursiveReflection(reflection, reflectionLevel + 1);
            fromThisLevel.r = (1.0 - ks.r) * fromLight.r + ks.r * fromReflection.r;
            fromThisLevel.g = (1.0 - ks.g) * fromLight.g + ks.g * fromReflection.g;
            fromThisLevel.b = (1.0 - ks.b) * fromLight.b + ks.b * fromReflection.b;

            return fromThisLevel;
        }
    }
    
}

/*
 Get color from each ray
 */
Color rayTraceColor(int x, int y)
{
    Color color;
    
    rayTracer ray = fireRays(x, y);
    
    // recursive reflection from 0 level
    color = recursiveReflection(ray, 0);
    
    return color;
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

//MODIFY THIS FUNCTION
void draw_scene()
{
  unsigned int x,y;
  //simple output
  for(x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(y=0;y < HEIGHT;y++)
    {
        // plot_pixel(x,y,x%256,y%256,(x+y)%256);
        Color rayColor = rayTraceColor(x, y);
        rayColor.clamp();
        // add the global ambient color once
        Color temp = rayColor.addLight(ambient_light);
        rayColor = temp;
        rayColor.clamp();
        plot_pixel(x, y, rayColor.r * 255, rayColor.g * 255, rayColor.b * 255);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

void parse_check(char *expected,char *found)
{
  if(strcasecmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(strcasecmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(strcasecmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(strcasecmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void display()
{

}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
      draw_scene();
      if(mode == MODE_JPEG)
	save_jpg();
    }
  once=1;
}

void keybutton(unsigned char key, int x, int y)
{
    switch (key){
        // 's' key is used to control: soft shadows.
        case 's':
            g_SoftShadows = !g_SoftShadows;
            if(g_SoftShadows){
                printf ("Generating scene with soft shadows...");
            }else{
                printf ("Generating scene with hard shadows...");
            }
            break;
            
        // 'q' 'w' 'e' 'r' keys are used to control: number of reflection times.
        case 'q':
            numReflection = 0;
            printf ("Generating scene only with primary ray...");
            break;
        case 'w':
            numReflection = 1;
            printf ("Generating scene with 1 times reflection...");
            break;
        case 'e':
            numReflection = 2;
            printf ("Generating scene with 2 times reflections...");
            break;
        case 'r':
            numReflection = 3;
            printf ("Generating scene with 3 times reflections...");
            break;
        
    default:
         break;
    }
    
    // refresh scene and save jpg images
    draw_scene();
    if(mode == MODE_JPEG) save_jpg();
}


int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      filename = argv[2];
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
    
  /* Due to issues of keyboard mapping on MAC, I use glutKeyboardFunc for g_ControlState switch */
  glutKeyboardFunc(keybutton);
    
  init();
  glutMainLoop();
}
