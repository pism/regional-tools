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

  int Mx() const {
    return m_Mx;
  }

  int My() const {
    return m_My;
  }

  double spacing() const {
    return m_spacing;
  }

  double dx() const {
    return m_dx;
  }

  double dy() const {
    return m_dy;
  }

  double x(size_t n) const {
    return m_x[n];
  }

  double y(size_t n) const {
    return m_y[n];
  }

private:
  double *m_x, *m_y, *m_z;
  double m_spacing, m_dx, m_dy;
  int m_Mx, m_My;

  double m_one_over_dx, m_one_over_dy;
  void get_corner_values(int i, int j,
                         double &A, double &B, double &C, double &D);
  Array2D<double> m_elevation;
};

inline int DEM::find_cell(const double *position, int &i, int &j) {
  i = floor((position[0] - m_x[0]) * m_one_over_dx);
  j = floor((position[1] - m_y[0]) * m_one_over_dy);

  // bail if we ended up outside the grid
  if (i < 0 || i + 1 > m_Mx - 1 || j < 0 || j + 1 > m_My - 1) {
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
  A = m_elevation(i,     j);
  B = m_elevation(i,     j + 1);
  C = m_elevation(i + 1, j + 1);
  D = m_elevation(i + 1, j);
}

#endif /* _DEM_H_ */
