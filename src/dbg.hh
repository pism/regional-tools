#ifndef _DBG_H_
#define _DBG_H_

enum MASK_VALUES {NO_VALUE = -2, ICE_FREE = -1};

int initialize_mask(int Mx, int My, double *thickness, int* mask);

int upslope_area(double *x, int Mx, double *y, int My, double *z, int *mask, bool output);

int accumulated_flow(double *x, int Mx, double *y, int My, double *z, double *my_mask, int n_samples);

#endif /* _DBG_H_ */
