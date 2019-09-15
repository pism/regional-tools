#include "dbg_internal.hh"
#include "DEM.hh"
#include <map>

/*!
 * Launch a particle from the starting point `i_start`, `j_start` and
 * track its trajectory as is advected by the velocity field `v =
 * -grad(F)`, where `F` is provided via `ctx`.
 *
 *
 */
static int streamline(dbg_context ctx, int i_start, int j_start,
                      Array2D<int> &old_mask, Array2D<int> &new_mask) {
  DEM *dem = (DEM*)ctx.system.params;

  int mask_counter = 0,
    n_max = (dem->Mx + dem->My) * ctx.steps_per_cell,
    i = i_start, j = j_start,
    i_old, j_old,
    status;

  double step_length = dem->spacing / ctx.steps_per_cell,
    position[2],
    err[2],
    gradient[2],
    elevation,
    gradient_magnitude;

  std::map<int,int> values;

  int mask_value = old_mask(i_start, j_start);

  // stop if the current cell already has a value assigned
  if (mask_value > 0 || mask_value == ICE_FREE)
    return 0;

  position[0] = dem->x[i_start] + dem->dx * 0.5;
  position[1] = dem->y[j_start] + dem->dy * 0.5;

  dem->evaluate(position, &elevation, NULL);

  // if there is no ice or we're below the minimum elevation, we're done
  if (elevation < ctx.min_elevation)
    return 0;

  // if we're above the maximum elevation, wait.
  if (elevation > ctx.max_elevation)
    return 1;

  for (int step_counter = 0; step_counter < n_max; ++step_counter) {

    i_old = i; j_old = j;
    status = dem->find_cell(position, i, j);

    if (status != 0)
      break;

    mask_value = old_mask(i, j);

    if (mask_value == ICE_FREE)   // ice-free
      break;

    if ((i != i_old || j != j_old) && (mask_value > 0)) {
      values[mask_value]++;
      mask_counter++;

      if (mask_counter == ctx.path_length)
        break;
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

  // Find the mask value that appears more often than others.
  int most_frequent_mask_value = NO_VALUE;
  int number_of_occurences = 0;

  std::map<int,int>::iterator k;
  for (k = values.begin(); k != values.end(); ++k) {
    if (k->second > number_of_occurences) {
      number_of_occurences = k->second;
      most_frequent_mask_value = k->first;
    }
  }

  new_mask(i_start, j_start) = most_frequent_mask_value;

  if (most_frequent_mask_value == NO_VALUE)
    return 1;

  return 0;
}

int upslope_area(double *x, int Mx, double *y, int My, double *z, int *mask, bool output) {
  int remaining = 0;
  const double elevation_step = 10;

  DEM dem(x, Mx, y, My, z);

  Array2D<int> my_mask(Mx, My), new_mask(Mx, My);
  my_mask.wrap(mask);
  if (new_mask.allocate() != 0)
    return -1;

#pragma omp parallel default(shared)
  {
    gsl_odeiv2_system system = {right_hand_side, NULL, 2, &dem};
    gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, 2);
    dbg_context ctx = {system, step, 2, 5, 0, elevation_step};

#pragma omp for schedule(dynamic)
    for (int j = 0; j < My; j++)
      for (int i = 0; i < Mx; i++)
        new_mask(i, j) = my_mask(i, j);

    do {
#pragma omp for schedule(dynamic)
      for (int j = 0; j < My; j++)
        for (int i = 0; i < Mx; i++)
          my_mask(i, j) = new_mask(i, j);

#pragma omp single
      remaining = 0;

#pragma omp for schedule(dynamic) reduction(+:remaining)
      for (int j = 0; j < My; j++) { // Note: traverse in the optimal order
        for (int i = 0; i < Mx; i++) {
          remaining += streamline(ctx, i, j, my_mask, new_mask);
        }
      }

      ctx.min_elevation = ctx.max_elevation;
      ctx.max_elevation += elevation_step;

    } while (remaining > 0);

    gsl_odeiv2_step_free(step);

#pragma omp for schedule(dynamic)
    for (int j = 0; j < My; j++)
      for (int i = 0; i < Mx; i++)
	my_mask(i, j) = new_mask(i, j);

  } // end of the parallel block

  return 0;
}
