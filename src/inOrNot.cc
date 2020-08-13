#include <iostream>
#include <stdio.h>
/* Determine whether the point (x, y) is in a polygon with vertices of
 * ploy_sides */
/* para: poly_sides    the number of vertices of the polygon
**    poly_x    the X-axis coordinates of each vertex of the polygon
**    poly_y     the Y-axis coordinates of each vertex of the polygon
**    x    X coordinate of test point
**    Y    Y-axis coordinate of test point*/
/* Return value: Return 0 means not inside the polygon, return 1 means inside
 * the polygon */
/* 在多边形各边上的点默认不在多边形内部 */

int inOrNot(int poly_sides, double* poly_X, double* poly_Y, double x, double y) {
  int i, j;
  j = poly_sides - 1;
  int res = 0;
  for (i = 0; i < poly_sides; i++) {
    //对每一条边进行遍历，该边的两个端点，有一个必须在待检测点(x,y)的左边，且两个点中，有一个点的y左边比p.y小，另一个点的y比p.y大。
    if ((poly_Y[i] < y && poly_Y[j] >= y || poly_Y[j] < y && poly_Y[i] >= y) &&
        (poly_X[i] <= x || poly_X[j] <= x)) {
      //用水平的直线与该边相交，求交点的x坐标。
      res ^= ((poly_X[i] + (y - poly_Y[i]) / (poly_Y[j] - poly_Y[i]) *
                               (poly_X[j] - poly_X[i])) < x);
    }
    j = i;
  }
  return res;
}