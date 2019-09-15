#include "DEM.hh"
#include <cmath>

DEM::DEM(double *x, int Mx, double *y, int My, double *z)
  : m_elevation(Mx, My, z) {
  m_x  = x;
  m_Mx = Mx;

  m_y  = y;
  m_My = My;

  m_z = z;

  // We assume that the grid is uniform.
  m_dx = x[1] - x[0];
  m_dy = y[1] - y[0];
  m_spacing = m_dx > m_dy ? m_dx : m_dy;

  m_one_over_dx = 1.0 / m_dx;
  m_one_over_dy = 1.0 / m_dy;
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
    delta_x = position[0] - m_x[i],
    delta_y = position[1] - m_y[j];

  // surface elevation
  if (elevation != NULL) {
    double
      alpha = m_one_over_dx * delta_x,
      beta  = m_one_over_dy * delta_y;

    *elevation = ( (1 - alpha) * (1 - beta) * A +
                   (1 - alpha) *      beta  * B +
                   alpha       *      beta  * C +
                   alpha       * (1 - beta) * D );
  }

  // the gradient
  if (gradient != NULL) {
    double gamma = m_one_over_dx * m_one_over_dy * (A + C - B - D);

    gradient[0] = (D - A) * m_one_over_dx + delta_x * gamma;
    gradient[1] = (B - A) * m_one_over_dy + delta_y * gamma;
  }

} // end of DEM::evaluate()
