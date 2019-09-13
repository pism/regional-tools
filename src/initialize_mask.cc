#include <string.h>             // memset
#include "Array2D.hh"
#include "dbg.hh"               // MASK_VALUES

int initialize_mask(int Mx, int My, double *thickness, int* output) {
  Array2D<double> THK(Mx, My);
  Array2D<int> mask(Mx, My), tmp(Mx, My);

  THK.wrap(thickness);
  mask.wrap(output);

  if (tmp.allocate() != 0)
    return 1;

  memset(mask.data(), 0, Mx*My*sizeof(int));

  double thk_eps = 1;

  int marker = 1;

  for (int i = 0; i < Mx; ++i) {
    for (int j = 0; j < My; ++j) {

      if (i == 0 || i == Mx - 1 ||
          j == 0 || j == My - 1) {
        tmp(i, j) = ICE_FREE;
        continue;
      }

      double thk = THK(i,     j),
        thk_w    = THK(i - 1, j),
        thk_nw   = THK(i - 1, j + 1),
        thk_n    = THK(i,     j + 1),
        thk_ne   = THK(i + 1, j + 1),
        thk_e    = THK(i + 1, j),
        thk_se   = THK(i + 1, j - 1),
        thk_s    = THK(i,     j - 1),
        thk_sw   = THK(i - 1, j - 1);

      if (thk > thk_eps) {
        // icy cell

        if (thk_w <= thk_eps || thk_nw <= thk_eps || thk_n <= thk_eps || thk_ne <= thk_eps ||
            thk_e <= thk_eps || thk_se <= thk_eps || thk_s <= thk_eps || thk_sw <= thk_eps) {
          // ice margin
          tmp(i, j) = ++marker;
        } else {
          // interior ice
          tmp(i, j) = NO_VALUE;
        }

      } else {
        // ice-free
        tmp(i, j) = ICE_FREE;
      }

    } // inner for loop
  } // outer for loop

  // second pass
  for (int i = 0; i < Mx; ++i) {
    for (int j = 0; j < My; ++j) {

      if (i == 0 || i == Mx - 1 ||
          j == 0 || j == My - 1) {
        mask(i, j) = tmp(i, j);
        continue;
      }

      double m  = tmp(i,     j),
        mask_w  = tmp(i - 1, j),
        mask_nw = tmp(i - 1, j + 1),
        mask_n  = tmp(i,     j + 1),
        mask_ne = tmp(i + 1, j + 1),
        mask_e  = tmp(i + 1, j),
        mask_se = tmp(i + 1, j - 1),
        mask_s  = tmp(i,     j - 1),
        mask_sw = tmp(i - 1, j - 1);

      double thk = THK(i,     j),
        thk_w    = THK(i - 1, j),
        thk_nw   = THK(i - 1, j + 1),
        thk_n    = THK(i,     j + 1),
        thk_ne   = THK(i + 1, j + 1),
        thk_e    = THK(i + 1, j),
        thk_se   = THK(i + 1, j - 1),
        thk_s    = THK(i,     j - 1),
        thk_sw   = THK(i - 1, j - 1);

      // ice-free cell next to an icy cell
      if (thk < thk_eps &&
          (thk_w >= thk_eps || thk_nw >= thk_eps || thk_n >= thk_eps || thk_ne >= thk_eps ||
           thk_e >= thk_eps || thk_se >= thk_eps || thk_s >= thk_eps || thk_sw >= thk_eps)) {

        if (mask_w > 0)
          mask(i, j) =  mask_w;
        else if (mask_nw > 0)
          mask(i, j) =  mask_nw;
        else if (mask_n > 0)
          mask(i, j) =  mask_n;
        else if (mask_ne > 0)
          mask(i, j) =  mask_ne;
        else if (mask_e > 0)
          mask(i, j) =  mask_e;
        else if (mask_se > 0)
          mask(i, j) =  mask_se;
        else if (mask_s > 0)
          mask(i, j) =  mask_s;
        else if (mask_sw > 0)
          mask(i, j) =  mask_sw;
      } else {
        mask(i, j) =  m;
      }

    } // inner for loop
  } // outer for loop

  return 0;
}
