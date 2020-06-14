#ifndef DRAW_H
#define DRAW_H

#include "matrix.h"
#include "ml6.h"
#include "symtab.h"

// Shading Type Constants
#define WIREFRAME 0
#define FLAT 1
#define GOURAUD 2
#define PHONG 3

// Scanline
void draw_scanline( double x0, double z0, double x1, double z1, int y, double offx,
                    color c0, color c1, screen s, zbuffer zb, 
                    double * view, double light[2][3], color ambient,
                    struct constants * reflect,
                    double * n0, double * n1, int type);
void scanline_convert(  struct matrix * points, int col, 
                        screen s, zbuffer zbuff, 
                        double * view, double light[2][3], color ambient,
                        struct constants * reflect, 
                        double norms[3][3], int type);

// Polygon organization
void add_polygons( struct matrix * polygons,
                   double x0, double y0, double z0,
                   double x1, double y1, double z1,
                   double x2, double y2, double z2);
void draw_polygons( struct matrix * polygons, screen s, zbuffer zb, 
                    double * view, double light[2][3], color ambient,
                    struct constants * reflect, int type, struct matrix * vns);

// Advanced shapes
// 3D shapes
void add_mesh(struct matrix * polygons, struct matrix * vns, char * filename);
void add_box( struct matrix * edges,
              double x, double y, double z,
              double width, double height, double depth );
void add_sphere( struct matrix * edges,
                 double cx, double cy, double cz,
                 double r, int step );
struct matrix * generate_sphere(double cx, double cy, double cz,
                                double r, int step );
void add_torus( struct matrix * edges,
                double cx, double cy, double cz,
                double r1, double r2, int step );
struct matrix * generate_torus( double cx, double cy, double cz,
                                double r1, double r2, int step );

// 2D Curves
void add_circle(struct matrix * points,
                double cx, double cy, double cz,
                double r, int step );
void add_curve( struct matrix *points,
                double x0, double y0,
                double x1, double y1,
                double x2, double y2,
                double x3, double y3,
                int step, int type );

void add_point(struct matrix * points, double x, double y, double z);
void add_edge(struct matrix * points,
	       double x0, double y0, double z0,
	       double x1, double y1, double z1);

void draw_lines(struct matrix * points, screen s, zbuffer zb, color c);
void draw_line( int x0, int y0, double z0, int x1, int y1, double z1, 
                screen s, zbuffer zb, color c);

// My funcs
void swap (double *a, double *b);
void change_color(color * c, int r, int g, int b);

#endif
