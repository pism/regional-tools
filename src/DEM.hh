#include <cmath>

#ifndef _DEM_H_
#define _DEM_H_

#include "Array2D.hh"

// A DEM on a regular cartesian grid. Uses bilinear interpolation.
class DEM {
public:
  DEM(double *x, int Mx, double *y, int My, double *z);
  ~DEM() {}

  int find_cell(const double *position, int &i, int &j);
  void evaluate(const double *position, double *elevation, double *f);

  double *x, *y, *z;
  double spacing, dx, dy;
  int Mx, My;
protected:
  double one_over_dx, one_over_dy;
  void get_corner_values(int i, int j,
                         double &A, double &B, double &C, double &D);
  Array2D<double> elevation;
};

inline int DEM::find_cell(const double *position, int &i, int &j) {
  i = floor((position[0] - x[0]) * one_over_dx);
  j = floor((position[1] - y[0]) * one_over_dy);

  // bail if we ended up outside the grid
  if (i < 0 || i + 1 > Mx - 1 || j < 0 || j + 1 > My - 1) {
    i = j = -1;
    return 1;
  }

  return 0;
}

inline void DEM::get_corner_values(int i, int j,
                                   double &A, double &B, double &C, double &D) {
  // Get the surface elevation at grid corners (arranged like so):
  //
  //   ^ y
  //   |
  //   |
  //   B-----C
  //   |     |
  //   | *   |   x
  // --A-----D---->
  //   |
  A = elevation(i,     j);
  B = elevation(i,     j + 1);
  C = elevation(i + 1, j + 1);
  D = elevation(i + 1, j);
}

#endif /* _DEM_H_ */
