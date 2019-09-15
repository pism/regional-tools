#include "dbg_internal.hh"
#include "DEM.hh"

static int streamline2(dbg_context ctx, int i_start, int j_start, Array2D<double> mask,
                       int n_samples) {
  DEM *dem = (DEM*)ctx.system.params;

  int n_max = (dem->Mx + dem->My) * ctx.steps_per_cell,
    i = i_start, j = j_start,
    i_old, j_old,
    status;

  double step_length = dem->spacing / ctx.steps_per_cell,
    position[2],
    err[2],
    gradient[2],
    gradient_magnitude;

  // stop if the current cell already has a value assigned
  if (mask(i_start, j_start) == ICE_FREE)
    return 0;

  // try M*M points inside the cell as starting points
  int M = n_samples;
  double step_x = dem->dx / (M + 1),
    step_y = dem->dy / (M + 1);

  for (int m = 1; m < M + 1; ++m) {
    for (int n = 1; n < M + 1; ++n) {

      position[0] = dem->x[i_start] + step_x * m;
      position[1] = dem->y[j_start] + step_y * n;

      for (int step_counter = 0; step_counter < n_max; ++step_counter) {

        i_old = i; j_old = j;
        status = dem->find_cell(position, i, j);

        if (status != 0)
          break;

        if (mask(i, j) == ICE_FREE)   // ice-free
          break;

        if (i != i_old || j != j_old) {
#pragma omp atomic
          mask(i, j) += 1.0 / (M * M);
        }

        dem->evaluate(position, NULL, gradient);

        gradient_magnitude = sqrt(gradient[0]*gradient[0] + gradient[1]*gradient[1]);

        // take a step
        status = gsl_odeiv2_step_apply(ctx.step,
                                       0,         // starting time (irrelevant)
                                       step_length / gradient_magnitude, // step size (units of time)
                                       position, err, NULL, NULL, &ctx.system);

        if (status != GSL_SUCCESS) {
          printf ("error, return value=%d\n", status);
          break;
        }

      } // time-stepping loop (step_counter)

    }
  }
  return 0;
}

int accumulated_flow(double *x, int Mx, double *y, int My, double *z, double *my_mask, int n_samples) {
  DEM dem(x, Mx, y, My, z);

  Array2D<double> mask(Mx, My);
  mask.wrap(my_mask);

#pragma omp parallel default(shared)
  {
    gsl_odeiv2_system system = {right_hand_side, NULL, 2, &dem};
    gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, 2);
    dbg_context ctx = {system, step, 2, // steps per cell
                       0, 0, 0};

#pragma omp for schedule(dynamic)
    for (int j = 0; j < My; j++) { // Note: traverse in the optimal order
      for (int i = 0; i < Mx; i++) {
        streamline2(ctx, i, j, mask, n_samples);
      }
    }

    gsl_odeiv2_step_free(step);

  } // end of the parallel block

  return 0;
}


