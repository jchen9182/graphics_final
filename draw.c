#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"
#include "gmath.h"
#include "symtab.h"

// scanline helper functions
color calc_color(color a, double * n, int counter) {
    // Minimize rounding error
    color c;
    double red = counter * n[0];
    int r = (red - floor(red) > 0.5) ? ceil(red) : floor(red);
    double green = counter * n[1];
    int g = (green - floor(green) > 0.5) ? ceil(green) : floor(green);
    double blue = counter * n[2];
    int b = (blue - floor(blue) > 0.5) ? ceil(blue) : floor(blue);

    c.red = a.red + r;
    c.green = a.green + g;
    c.blue = a.blue + b;
    return c;
}
void swapc(color *a, color *b) {
    color temp = *a;
    *a = *b;
    *b = temp;
}
void print_color(color c) { // debugging help
    printf("%u %u %u\n", c.red, c.green, c.blue);
}

void copy_array(double * a, double * b) { // copies b to a
    a[0] = b[0];
    a[1] = b[1];
    a[2] = b[2];
}
void swap_norm(double * a, double * b) {
    double temp[3];
    copy_array(temp, a);
    copy_array(a, b);
    copy_array(b, temp);
}
void add_norm(double * a, double * b) {
    a[0] += b[0];
    a[1] += b[1];
    a[2] += b[2];
}
void print_norm(double * n) {
    printf("%f %f %f\n", n[0], n[1], n[2]);
}
/*======== void draw_scanline() ==========
  Inputs: struct matrix *points
          int i
          screen s
          zbuffer zb
          color c
  Line algorithm specifically for horizontal scanlines
  ====================*/
void draw_scanline( double x0, double z0, double x1, double z1, int y, double offx,
                    color c0, color c1, screen s, zbuffer zb, 
                    double * view, double light[2][3], color ambient,
                    struct constants * reflect,
                    double * n0, double * n1, int type) {
    if (x0 > x1) {
        swap(&x0, &x1);
        swap(&z0, &z1);
        if (type == GOURAUD) swapc(&c0, &c1);
        if (type == PHONG) swap_norm(n0, n1);
    }

    int x = ceil(x0);
    double dist = x1 - x0 + 1;
    
    double mz = dist > 0 ? (z1 - z0) / dist : 0;
    double z = z0 + mz * offx;

    double mc[3];
    if (type == GOURAUD) {
        mc[0] = dist > 0 ? (c1.red - c0.red) / dist : 0;
        mc[1] = dist > 0 ? (c1.green - c0.green) / dist : 0;
        mc[2] = dist > 0 ? (c1.blue - c0.blue) / dist : 0;
    }

    double mn[3];
    if (type == PHONG) {
        mn[0] = dist > 0 ? (n1[0] - n0[0]) / dist : 0;
        mn[1] = dist > 0 ? (n1[1] - n0[1]) / dist : 0;
        mn[2] = dist > 0 ? (n1[2] - n0[2]) / dist : 0;
    }

    color c = c0;
    while (x < ceil(x1)) {
        if (type == PHONG) c = get_lighting(n0, view, ambient, light, reflect);

        plot(s, zb, c, x, y, z);

        z += mz;
        x++;

        if (type == GOURAUD) {
            c = calc_color(c0, mc, x - ceil(x0));
        }
        else if (type == PHONG) {
            add_norm(n0, mn);
        }
    }
}

// scanline_convert helper functions
void swap(double *a, double *b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}
/*======== void scanline_convert() ==========
  Inputs: struct matrix *points
          int i
          screen s
          zbuffer zb
  Returns:
  Fills in polygon i by drawing consecutive horizontal (or vertical) lines.
  Color should be set differently for each polygon.
  Includes Greg's pixel perfect scanning
  ====================*/
void scanline_convert(  struct matrix * points, int col, 
                        screen s, zbuffer zbuff, 
                        double * view, double light[2][3], color ambient,
                        struct constants * reflect, 
                        double norms[3][3], int type) {
    double ** matrix = points -> m;
    double xb = matrix[0][col];
    double xm = matrix[0][col + 1];
    double xt = matrix[0][col + 2];
    double yb = matrix[1][col];
    double ym = matrix[1][col + 1];
    double yt = matrix[1][col + 2];
    double zb = matrix[2][col];
    double zm = matrix[2][col + 1];
    double zt = matrix[2][col + 2];

    double nb[3], nm[3], nt[3];
    if (type > 1) {
        copy_array(nb, norms[0]);
        copy_array(nm, norms[1]);
        copy_array(nt, norms[2]);
    }

    if (yb > ym) {
        swap(&xb, &xm);
        swap(&yb, &ym);
        swap(&zb, &zm);
        if (type > 1) swap_norm(nb, nm);
    }
    if (ym > yt) {
        swap(&xm, &xt);
        swap(&ym, &yt);
        swap(&zm, &zt);
        if (type > 1) swap_norm(nm, nt);
    }
    if (yb > ym) {
        swap(&xb, &xm);
        swap(&yb, &ym);
        swap(&zb, &zm);
        if (type > 1) swap_norm(nb, nm);
    }

    color cb, cm, ct;
    if (type == GOURAUD) {
        cb = get_lighting(nb, view, ambient, light, reflect);
        cm = get_lighting(nm, view, ambient, light, reflect);
        ct = get_lighting(nt, view, ambient, light, reflect);
    }

    double dist0 = yt - yb + 1;
    double dist1 = ym - yb + 1;
    double dist2 = yt - ym + 1;

    double mx0 = dist0 > 0 ? (xt - xb) / dist0 : 0;
    double mx1 = dist1 > 0 ? (xm - xb) / dist1 : 0;
    double mx2 = dist2 > 0 ? (xt - xm) / dist2 : 0;
    double mz0 = dist0 > 0 ? (zt - zb) / dist0 : 0;
    double mz1 = dist1 > 0 ? (zm - zb) / dist1 : 0;
    double mz2 = dist2 > 0 ? (zt - zm) / dist2 : 0;

    double mc0[3];
    double mc1[3];
    double mc2[3];
    if (type == GOURAUD) {
        mc0[0] = dist0 > 0 ? (ct.red - cb.red) / dist0 : 0;
        mc0[1] = dist0 > 0 ? (ct.green - cb.green) / dist0 : 0;
        mc0[2] = dist0 > 0 ? (ct.blue - cb.blue) / dist0 : 0;
        mc1[0] = dist1 > 0 ? (cm.red - cb.red) / dist1 : 0;
        mc1[1] = dist1 > 0 ? (cm.green - cb.green) / dist1 : 0;
        mc1[2] = dist1 > 0 ? (cm.blue - cb.blue) / dist1 : 0;
        mc2[0] = dist2 > 0 ? (ct.red - cm.red) / dist2 : 0;
        mc2[1] = dist2 > 0 ? (ct.green - cm.green) / dist2 : 0;
        mc2[2] = dist2 > 0 ? (ct.blue - cm.blue) / dist2 : 0;
    }

    double mn0[3];
    double mn1[3];
    double mn2[3];
    if (type == PHONG) {
        mn0[0] = dist0 > 0 ? (nt[0] - nb[0]) / dist0 : 0;
        mn0[1] = dist0 > 0 ? (nt[1] - nb[1]) / dist0 : 0;
        mn0[2] = dist0 > 0 ? (nt[2] - nb[2]) / dist0 : 0;
        mn1[0] = dist1 > 0 ? (nm[0] - nb[0]) / dist1 : 0;
        mn1[1] = dist1 > 0 ? (nm[1] - nb[1]) / dist1 : 0;
        mn1[2] = dist1 > 0 ? (nm[2] - nb[2]) / dist1 : 0;
        mn2[0] = dist2 > 0 ? (nt[0] - nm[0]) / dist2 : 0;
        mn2[1] = dist2 > 0 ? (nt[1] - nm[1]) / dist2 : 0;
        mn2[2] = dist2 > 0 ? (nt[2] - nm[2]) / dist2 : 0;
    }

    double offy0 = ceil(yb) - yb;
    double offy1 = ceil(ym) - ym;

    double x0 = xb + mx0 * offy0;
    double x1 = xb + mx1 * offy0;
    double x2 = xm + mx2 * offy1;
    double z0 = zb + mz0 * offy0;
    double z1 = zb + mz1 * offy0;
    double z2 = zm + mz2 * offy1;
    int y = ceil(yb);
    
    color c0, c1, cinit;
    if (type == GOURAUD) {
        c0 = cb;
        c1 = cb;
        cinit = cb;
    }

    double n0[3], n1[3];
    if (type == PHONG) {
        copy_array(n0, nb);
        copy_array(n1, nb);
    }

    int toggle = 1;
    while (y < ceil(yt)) {
        double offx;

        if (y == ceil(ym) && toggle) {
            x1 = x2;
            z1 = z2;
            mx1 = mx2;
            mz1 = mz2;

            if (type == GOURAUD) {
                c1 = cm;
                cinit = cm;
                copy_array(mc1, mc2);
            }

            if (type == PHONG) {
                copy_array(n1, nm);
                copy_array(mn1, mn2);
            }

            toggle = 0;
        }

        if (x0 > x1) {
            offx = ceil(x1) - x1;
        }
        else offx = ceil(x0) - x0;

        // dummy normals
        double dn0[3] = {0, 0, 0};
        double dn1[3] = {0, 0, 0};

        if (type == FLAT) {
            color c; // dummy color
            color cflat = get_lighting(norms[0], view, ambient, light, reflect);
            draw_scanline(x0, z0, x1, z1, y, offx, cflat, c, s, zbuff, 
            view, light, ambient, reflect, dn0, dn1, type);
        }
        else if (type == GOURAUD) {
            draw_scanline(x0, z0, x1, z1, y, offx, c0, c1, s, zbuff, 
            view, light, ambient, reflect, dn0, dn1, type);
        }
        else { // type = PHONG
            draw_scanline(x0, z0, x1, z1, y, offx, c0, c1, s, zbuff, 
            view, light, ambient, reflect, n0, n1, type);
        }

        x0 += mx0;
        x1 += mx1;
        z0 += mz0;
        z1 += mz1;
        y++;

        if (type == GOURAUD) {
            int counter0 = y - ceil(yb);
            int counter1 = (y >= ceil(ym)) ? y - ceil(ym) : counter0;

            c0 = calc_color(cb, mc0, counter0);
            c1 = calc_color(cinit, mc1, counter1);
        }
        else if (type == PHONG) {
            add_norm(n0, mn0);
            add_norm(n1, mn1);
        }
    }
}
/*======== void add_polygon() ==========
  Inputs:   struct matrix *polygons
            double x0
            double y0
            double z0
            double x1
            double y1
            double z1
            double x2
            double y2
            double z2
  Returns:
  Adds the vertices (x0, y0, z0), (x1, y1, z1)
  and (x2, y2, z2) to the polygon matrix. They
  define a single triangle surface.
  ====================*/
void add_polygon(struct matrix * polygons, 
                 double x0, double y0, double z0, 
                 double x1, double y1, double z1, 
                 double x2, double y2, double z2) {
    add_point(polygons, x0, y0, z0);
    add_point(polygons, x1, y1, z1);
    add_point(polygons, x2, y2, z2);
}

// draw_polygons helper function
int compare(double ** matrix, int thisCol, int compareCol) {
    double x = matrix[0][thisCol];
    double y = matrix[1][thisCol];
    double z = matrix[2][thisCol];

    for (int i = compareCol; i < compareCol + 3; i++) {
        double thisX = matrix[0][i];
        double thisY = matrix[1][i];
        double thisZ = matrix[2][i];

        if (x == thisX && y == thisY && z == thisZ) return 1;
    }
    return 0;
}
/*======== void draw_polygons() ==========
  Inputs:   struct matrix *polygons
            screen s
            color c
  Returns:
  Goes through polygons 3 points at a time, drawing
  lines connecting each points to create bounding triangles
  ====================*/
void draw_polygons( struct matrix * polygons, screen s, zbuffer zb, 
                    double * view, double light[2][3], color ambient,
                    struct constants * reflect, int type, struct matrix * vns) {
    int lastcol = polygons -> lastcol;
    double ** matrix = polygons -> m;

    if (lastcol < 3) {
        printf("Need at least 3 points to draw a polygon!\n");
        return;
    }

    if (type < 2) {
        int counter = 0;
        for (int col = 0; col < lastcol - 2; col += 3) {
            double * normal = calculate_normal(polygons, col);

            if (type == WIREFRAME) {
                double x0 = matrix[0][col];
                double y0 = matrix[1][col];
                double x1 = matrix[0][col + 1];
                double y1 = matrix[1][col + 1];
                double x2 = matrix[0][col + 2];
                double y2 = matrix[1][col + 2];
                color c;
                c.red = 0;
                c.green = 0;
                c.blue = 0;

                draw_line(x0, y0, 0, x1, y1, 0, s, zb, c);
                draw_line(x1, y1, 0, x2, y2, 0, s, zb, c);
                draw_line(x2, y2, 0, x0, y0, 0, s, zb, c);
            }
            else { // type = FLAT
                double norms[3][3] = {normal[0], normal[1], normal[2], 0, 0, 0, 0, 0, 0};
                if (normal[2] > 0) {
                    scanline_convert(polygons, col, s, zb, 
                    view, light, ambient, reflect,
                    norms, type);
                }
            }
        }
    }

    else {
        // Using predetermined normals only work in the original
        // orientation of an image. We have to calculate normals manually
        // otherwise. Change the > to a == for a rough fix
        if (vns -> lastcol > 0) { // mesh
            double norms[3][3];
            double ** normals = vns -> m;

            for (int col = 0; col < lastcol - 2; col += 3) {
                for (int i = 0; i < 3; i++) {
                    norms[i][0] = normals[0][col + i];
                    norms[i][1] = normals[1][col + i];
                    norms[i][2] = normals[2][col + i];
                }

                if (normals[2][col] > 0) {
                    scanline_convert(polygons, col, s, zb, 
                    view, light, ambient, reflect,
                    norms, type);
                }
            }
        }
        else { // regular shapes
            double norms[3][3];
            int counter = 0;

            for (int col = 0; col < lastcol; col++) {
                double averageNormal[3] = {0, 0, 0};

                for (int c = 0; c < lastcol - 2; c += 3) {
                    if (compare(matrix, col, c)) {
                        double * normal = calculate_normal(polygons, c);
                        if (normal[2] > 0) {
                            averageNormal[0] += normal[0];
                            averageNormal[1] += normal[1];
                            averageNormal[2] += normal[2];
                        }
                    }
                }
                normalize(averageNormal);
                copy_array(norms[counter], averageNormal);

                counter++;
                if (counter % 3 == 0 && col > 0) {
                    double * normal = calculate_normal(polygons, col - 2);

                    if (normal[2] > 0) {
                        scanline_convert(polygons, col - 2, s, zb, 
                        view, light, ambient, reflect,
                        norms, type);
                    }
                    counter = 0;
                }
            }
        }
    }
}

void add_mesh(struct matrix * polygons, struct matrix * vns, char * filename) {
    FILE * fp;
    char buffer[255];

    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Couldn't Open File.\n");
        return;
    }
    else {
        // Find number of v, vn
        int vnum = 0;
        int normnum = 0;
        while (fgets(buffer, 255, fp)) {
            if (buffer[0] == 'v') {
                if (buffer[1] == 'n') normnum++;
                else vnum++;
            }
        }

        // Store vertex and face data
        double vertices[vnum][3];
        double vertnorms[normnum][3];
        int vcounter = 0;
        int vncounter = 0;
        fseek(fp, 0, SEEK_SET);

        while (fgets(buffer, 255, fp)) {
            if (buffer[0] == 'v') {
                if (buffer[1] == 'n') {
                    double x, y, z;
                    
                    if (buffer[3] == ' ') { 
                        // Accounts for the extra space when there are negative numbers
                        sscanf(buffer, "vn  %le %le %le", &x, &y, &z);
                    }
                    else {
                        sscanf(buffer, "vn %le %le %le", &x, &y, &z);
                    }
                    vertnorms[vncounter][0] = x;
                    vertnorms[vncounter][1] = y;
                    vertnorms[vncounter][2] = z;

                    vncounter++;
                }
                else if (buffer[1] == ' '){
                    double x, y, z;

                    if (buffer[2] == ' ') {
                        sscanf(buffer, "v  %le %le %le", &x, &y, &z);
                    }
                    else {
                        sscanf(buffer, "v %le %le %le", &x, &y, &z);
                    }
                    vertices[vcounter][0] = x;
                    vertices[vcounter][1] = y;
                    vertices[vcounter][2] = z;

                    vcounter++;
                }
            }
            else if(buffer[0] == 'f') {
                // only works with v//vn and v/vt/vn formats
                int v0, vn0, vt0, v1, vn1, vt1, v2, vn2, vt2, v3, vn3, vt3;

                int slashCounter = 0;
                int hasTexture = 1;

                for (int i = 1; i < strlen(buffer); i++) {
                    if (buffer[i] == '/') slashCounter++;
                    if (buffer[i] == '/' && buffer[i] == buffer[i + 1])
                        hasTexture = 0;
                }

                if (slashCounter < 6) { // triangle faces
                    if (hasTexture) {
                        sscanf(buffer, "f %d/%d/%d %d/%d/%d %d/%d/%d", 
                        &v0, &vt0, &vn0, &v1, &vt1, &vn1, &v2, &vt2, &vn2);
                    }
                    else {
                        sscanf(buffer, "f %d//%d %d//%d %d//%d", 
                        &v0, &vn0, &v1, &vn1, &v2, &vn2);
                    }
                }
                else { // quadrilateral faces
                    if (hasTexture) {
                        sscanf(buffer, "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d", 
                        &v0, &vt0, &vn0, &v1, &vt1, &vn1, 
                        &v2, &vt2, &vn2, &v3, &vt3, &vn3);
                    }
                    else {
                        sscanf(buffer, "f %d//%d %d//%d %d//%d %d//%d", 
                        &v0, &vn0, &v1, &vn1, &v2, &vn2, &v3, &vn3);
                    }
                }

                // This only works for obj files where all the faces are at the end of the file
                add_polygon(polygons, 
                    vertices[v0 - 1][0], vertices[v0 - 1][1], vertices[v0 - 1][2],
                    vertices[v1 - 1][0], vertices[v1 - 1][1], vertices[v1 - 1][2],
                    vertices[v2 - 1][0], vertices[v2 - 1][1], vertices[v2 - 1][2]);
                add_polygon(vns, 
                    vertnorms[vn0 - 1][0], vertnorms[vn0 - 1][1], vertnorms[vn0 - 1][2],
                    vertnorms[vn1 - 1][0], vertnorms[vn1 - 1][1], vertnorms[vn1 - 1][2],
                    vertnorms[vn2 - 1][0], vertnorms[vn2 - 1][1], vertnorms[vn2 - 1][2]);
                if (slashCounter > 6) {
                    add_polygon(polygons, 
                    vertices[v0 - 1][0], vertices[v0 - 1][1], vertices[v0 - 1][2],
                    vertices[v2 - 1][0], vertices[v2 - 1][1], vertices[v2 - 1][2],
                    vertices[v3 - 1][0], vertices[v3 - 1][1], vertices[v3 - 1][2]);
                    add_polygon(vns, 
                    vertnorms[vn0 - 1][0], vertnorms[vn0 - 1][1], vertnorms[vn0 - 1][2],
                    vertnorms[vn2 - 1][0], vertnorms[vn2 - 1][1], vertnorms[vn2 - 1][2],
                    vertnorms[vn3 - 1][0], vertnorms[vn3 - 1][1], vertnorms[vn3 - 1][2]);
                }
            }
        }
    }
}

/*======== void add_box() ==========
  Inputs:   struct matrix * edges
            double x
            double y
            double z
            double width
            double height
            double depth
  add the points for a rectagular prism whose
  upper-left-front corner is (x, y, z) with width,
  height and depth dimensions.
  ====================*/
void add_box(struct matrix * polygons,
             double x, double y, double z,
             double width, double height, double depth) {
    double x1 = x + width;
    double y1 = y - height;
    double z1 = z - depth;

    // Left Face
    add_polygon(polygons, x, y, z1, x, y1, z1, x, y, z);
    add_polygon(polygons, x, y, z, x, y1, z1, x, y1, z);
    // // Right Face
    add_polygon(polygons, x1, y, z, x1, y1, z, x1, y, z1);
    add_polygon(polygons, x1, y, z1, x1, y1, z, x1, y1, z1);

    // Front Face
    add_polygon(polygons, x, y, z, x, y1, z, x1, y, z);
    add_polygon(polygons, x1, y, z, x, y1, z, x1, y1, z);
    // Back Face
    add_polygon(polygons, x1, y, z1, x1, y1, z1, x, y, z1);
    add_polygon(polygons, x, y, z1, x1, y1, z1, x, y1, z1);

    // Top Face
    add_polygon(polygons, x, y, z1, x, y, z, x1, y, z1);
    add_polygon(polygons, x1, y, z1, x, y, z, x1, y, z);
    // Bottom Face
    add_polygon(polygons, x, y1, z, x, y1, z1, x1, y1, z);
    add_polygon(polygons, x1, y1, z, x, y1, z1, x1, y1, z1);
}

/*======== void add_sphere() ==========
  Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r
            int step
  adds all the points for a sphere with center (cx, cy, cz)
  and radius r using step points per circle/semicircle.
  Since edges are drawn using 2 points, add each point twice,
  or add each point and then another point 1 pixel away.
  should call generate_sphere to create the necessary points
  ====================*/
void add_sphere(struct matrix * polygons,
                double cx, double cy, double cz,
                double r, int step) {
    struct matrix * sphere = generate_sphere(cx, cy, cz, r, step);
    double ** matrix = sphere -> m;
    
    for (int lat = 0; lat < step; lat++) {
        for (int longt = 0; longt < step; longt++) {
            int index = lat * (step + 1) + longt;

            int p0 = index;
            int p1 = index + 1;
            int p2 = (index + step) % (step * (step + 1));
            int p3 = (index + step + 1) % (step * (step + 1));

            double x0 = matrix[0][p0];
            double y0 = matrix[1][p0];
            double z0 = matrix[2][p0];

            double x1 = matrix[0][p1];
            double y1 = matrix[1][p1];
            double z1 = matrix[2][p1];

            double x2 = matrix[0][p2];
            double y2 = matrix[1][p2];
            double z2 = matrix[2][p2];

            double x3 = matrix[0][p3];
            double y3 = matrix[1][p3];
            double z3 = matrix[2][p3];

            add_polygon(polygons, x0, y0, z0, x3, y3, z3, x2, y2, z2);
            add_polygon(polygons, x0, y0, z0, x1, y1, z1, x3, y3, z3);
        }
    }
    free_matrix(sphere);
}

/*======== void generate_sphere() ==========
  Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r
            int step
  Returns: Generates all the points along the surface
           of a sphere with center (cx, cy, cz) and
           radius r using step points per circle/semicircle.
           Returns a matrix of those points
  ====================*/
struct matrix * generate_sphere(double cx, double cy, double cz,
                                double r, int step) {
    struct matrix * points = new_matrix(4, step * step);

    for (int p = 0; p < step; p++) {
        double phi = (2 * M_PI) * p / step;

        for (int t = 0; t <= step; t++) {
            double theta = M_PI * t / step;

            double x = r * cos(theta) + cx;
            double y = r * sin(theta) * cos(phi) + cy;
            double z = r * sin(theta) * sin(phi) + cz;
            add_point(points, x, y, z);
        }
    }

    return points;
}

/*======== void add_torus() ==========
  Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r1
            double r2
            double step
  Returns:
  adds all the points required for a torus with center (cx, cy, cz),
  circle radius r1 and torus radius r2 using step points per circle.
  should call generate_torus to create the necessary points
  ====================*/
void add_torus( struct matrix * polygons,
                double cx, double cy, double cz,
                double r1, double r2, int step) {
    struct matrix * torus = generate_torus(cx, cy, cz, r1, r2, step);
    double ** matrix = torus -> m;

    for (int lat = 0; lat < step; lat++) {
        for (int longt = 0; longt < step; longt++) {
            int index = lat * step + longt;

            int p0 = index;
            int p1 = index + 1;
            if (longt == step - 1) p1 = index - longt;
            int p2 = (index + step) % (step * step);
            int p3 = (p1 + step) % (step * step);

            double x0 = matrix[0][p0];
            double y0 = matrix[1][p0];
            double z0 = matrix[2][p0];

            double x1 = matrix[0][p1];
            double y1 = matrix[1][p1];
            double z1 = matrix[2][p1];

            double x2 = matrix[0][p2];
            double y2 = matrix[1][p2];
            double z2 = matrix[2][p2];

            double x3 = matrix[0][p3];
            double y3 = matrix[1][p3];
            double z3 = matrix[2][p3];

            add_polygon(polygons, x0, y0, z0, x2, y2, z2, x3, y3, z3);
            add_polygon(polygons, x0, y0, z0, x3, y3, z3, x1, y1, z1);
        }
    }
    free_matrix(torus);
}

/*======== void generate_torus() ==========
  Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r
            int step
  Returns: Generates all the points along the surface
           of a torus with center (cx, cy, cz),
           circle radius r1 and torus radius r2 using
           step points per circle.
           Returns a matrix of those points
  ====================*/
struct matrix * generate_torus( double cx, double cy, double cz,
                                double r1, double r2, int step) {
    struct matrix * points = new_matrix(4, step * step);
    double phi = 0;
    double theta = 0;

    for (int p = 0; p < step; p++) {
        phi = (2 * M_PI) * p / step;

        for (int t = 0; t < step; t++) {
            theta = (2 * M_PI) * t / step;

            double x = cos(phi) * (r1 * cos(theta) + r2) + cx;
            double y = r1 * sin(theta) + cy;
            double z = -1 * sin(phi) * (r1 * cos(theta) + r2) + cz;

            add_point(points, x, y, z);
        }
    }

    return points;
}

/*======== void add_circle() ==========
  Inputs:   struct matrix * edges
            double cx
            double cy
            double r
            double step
  Adds the circle at (cx, cy) with radius r to edges
  ====================*/
void add_circle(struct matrix *edges,
                double cx, double cy, double cz,
                double r, int step) {
    double angle = 0;
    double x0 = cx + r;
    double y0 = cy;

    for (int t = 0; t <= step; t++) {
        double x1 = r * cos(angle) + cx;
        double y1 = r * sin(angle) + cy;

        add_edge(edges, x0, y0, cz, x1, y1, cz);
        x0 = x1;
        y0 = y1;
        angle += (2 * M_PI) / step;
    }
}

/*======== void add_curve() ==========
Inputs:   struct matrix *edges
         double x0
         double y0
         double x1
         double y1
         double x2
         double y2
         double x3
         double y3
         double step
         int type
Adds the curve bounded by the 4 points passsed as parameters
of type specified in type (see matrix.h for curve type constants)
to the matrix edges
====================*/
void add_curve( struct matrix *edges,
                double x0, double y0,
                double x1, double y1,
                double x2, double y2,
                double x3, double y3,
                int step, int type) {
    struct matrix * xco = generate_curve_coefs(x0, x1, x2, x3, type);
    struct matrix * yco = generate_curve_coefs(y0, y1, y2, y3, type);
    double ** xm = xco -> m;
    double ** ym = yco -> m;
    double xold = x0;
    double yold = y0;

    for (int i = 0; i <= step; i++) {
        double xnew, ynew;
        double t = (double) i / step;

        xnew = t * (t * (xm[0][0] * t + xm[1][0]) + xm[2][0]) + xm[3][0];
        ynew = t * (t * (ym[0][0] * t + ym[1][0]) + ym[2][0]) + ym[3][0];

        add_edge(edges, xold, yold, 0, xnew, ynew, 0);
        xold = xnew;
        yold = ynew;
    }
}

/*======== void add_point() ==========
Inputs:   struct matrix * points
int x
int y
int z
Returns:
adds point (x, y, z) to points and increment points.lastcol
if points is full, should call grow on points
====================*/
void add_point(struct matrix * points, double x, double y, double z) {
    double ** matrix = points -> m;
    int cols = points -> cols;
    int lastcol = points -> lastcol;

    if (lastcol == cols) grow_matrix(points, cols + 100);

    matrix[0][lastcol] = x;
    matrix[1][lastcol] = y;
    matrix[2][lastcol] = z;
    matrix[3][lastcol] = 1;

    points -> lastcol++;
}

/*======== void add_edge() ==========
Inputs:   struct matrix * points
int x0, int y0, int z0, int x1, int y1, int z1
Returns:
add the line connecting (x0, y0, z0) to (x1, y1, z1) to points
should use add_point
====================*/
void add_edge(struct matrix * points,
	       double x0, double y0, double z0,
	       double x1, double y1, double z1) {
    add_point(points, x0, y0, z0);
    add_point(points, x1, y1, z1);
}

/*======== void draw_lines() ==========
Inputs:   struct matrix * points
screen s
color c
Returns:
Go through points 2 at a time and call draw_line to add that line
to the screen
====================*/
void draw_lines(struct matrix * points, screen s, zbuffer zb, color c) {
    int lastcol = points -> lastcol;

    if (lastcol < 2) {
        printf("Need at least 2 points to draw a line!\n");
        return;
    }

    for (int point = 0; point < lastcol - 1; point += 2) {
        int x0 = points -> m[0][point];
        int y0 = points -> m[1][point];
        double z0 = points -> m[2][point];
        int x1 = points -> m[0][point + 1];
        int y1 = points -> m[1][point + 1];
        double z1 = points -> m[2][point + 1];

        draw_line(x0, y0, z0, x1, y1, z1, s, zb, c);
    }
}

void draw_line( int x0, int y0, double z0, int x1, int y1, double z1, 
                screen s, zbuffer zb, color c) {
    int x, y, d, A, B;
    double z, mz;

    //swap points if going right -> left
    int xt, yt, zt;
    if (x0 > x1) {
        xt = x0;
        yt = y0;
        zt = z0;
        x0 = x1;
        y0 = y1;
        z0 = z1;
        x1 = xt;
        y1 = yt;
        z1 = zt;
    }

    x = x0;
    y = y0;
    z = z0;
    A = 2 * (y1 - y0);
    B = -2 * (x1 - x0);

    //octants 1 and 8
    if ( abs(x1 - x0) >= abs(y1 - y0) ) {
        mz = (z1 - z0) / (x1 - x0);
        //octant 1
        if ( A > 0 ) {
            d = A + B/2;

            while ( x < x1 ) {
                plot( s, zb, c, x, y, z );

                if ( d > 0 ) {
                    y+= 1;
                    d+= B;
                }

                x++;
                d+= A;
                z += mz;
            }

            plot( s, zb, c, x1, y1, z1 );
        }

        //octant 8
        else {
            d = A - B/2;

            while ( x < x1 ) {
                plot( s, zb, c, x, y, z );

                if ( d < 0 ) {
                    y-= 1;
                    d-= B;
                }

                x++;
                d+= A;
                z += mz;
            }

            plot( s, zb, c, x1, y1, z1 );
        }
    }

    //octants 2 and 7
    else {
        mz = (z1 - z0) / (y1 - y0);
        //octant 2
        if ( A > 0 ) {
            d = A/2 + B;

            while ( y < y1 ) {

                plot( s, zb, c, x, y, z );

                if ( d < 0 ) {
                    x+= 1;
                    d+= A;
                }

                y++;
                d+= B;
                z += mz;
            } 

            plot( s, zb, c, x1, y1, z1 );
        }

        //octant 7
        else {
            d = A/2 - B;

            while ( y > y1 ) {

                plot( s, zb, c, x, y, z );

                if ( d > 0 ) {
                    x+= 1;
                    d+= A;
                }

                y--;
                d-= B;
                z += mz;
            } 

            plot( s, zb, c, x1, y1, z1 );
        } 
    }
}