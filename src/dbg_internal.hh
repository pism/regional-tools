#ifndef _DBG_INTERNAL_H_
#define _DBG_INTERNAL_H_

// GSL stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "dbg.hh"
#include "DEM.hh"

struct dbg_context {
  gsl_odeiv2_system system;
  gsl_odeiv2_step *step;
  int steps_per_cell;
  // parameters below are used by the upslope area computation only
  int path_length;
  double min_elevation, max_elevation;
};

static int right_hand_side(double t, const double y[], double dydt[], void* params) {

  ((DEM*)params)->evaluate(y, NULL, dydt);

  dydt[0] *= -1;
  dydt[1] *= -1;

  return GSL_SUCCESS;

}

#endif /* _DBG_INTERNAL_H_ */
