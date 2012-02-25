#include "DEM.hh"
#include <cmath>

DEM::DEM(double *my_x, int my_Mx, double *my_y, int my_My, double *my_z)
  : elevation(my_Mx, my_My) {
  x  = my_x;
  Mx = my_Mx;

  y  = my_y;
  My = my_My;

  z   = my_z;

  // We assume that the grid is uniform.
  dx = x[1] - x[0];
  dy = y[1] - y[0];
  spacing = dx > dy ? dx : dy;

  one_over_dx = 1.0 / dx;
  one_over_dy = 1.0 / dy;

  elevation.wrap(z);
}

void DEM::evaluate(const double *position, double *elevation, double *gradient) {
  int ierr, i, j;
  double A, B, C, D;

  ierr = this->find_cell(position, i, j);

  // Pretend that outside the grid the surface is perfectly flat, the elevation
  // of the sea level (0) and ice-free.
  if (ierr != 0) {

    if (elevation != NULL)
      *elevation = 0;

    if (gradient != NULL)
      gradient[0] = gradient[1] = 0;

    return;
  }

  this->get_corner_values(i, j, A, B, C, D);

  double
    delta_x = position[0] - x[i],
    delta_y = position[1] - y[j];

  // surface elevation
  if (elevation != NULL) {
    double
      alpha = one_over_dx * delta_x,
      beta  = one_over_dy * delta_y;

    *elevation = ( (1 - alpha) * (1 - beta) * A +
                   (1 - alpha) *      beta  * B +
                   alpha       *      beta  * C +
                   alpha       * (1 - beta) * D );
  }

  // the gradient
  if (gradient != NULL) {
    double gamma = one_over_dx * one_over_dy * (A + C - B - D);

    gradient[0] = (D - A) * one_over_dx + delta_x * gamma;
    gradient[1] = (B - A) * one_over_dy + delta_y * gamma;
  }

} // end of DEM::evaluate()
