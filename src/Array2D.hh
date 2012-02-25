#include <cstdlib>

#ifndef _ARRAY2D_H_
#define _ARRAY2D_H_

/* A zero-cost 2D array class hiding storage order details. */
/* This is a way to have clean code and avoid using evil macros. */
template<class T>
class Array2D {
public:
  Array2D(int my_Mx, int my_My)
    : wrapper(false), private_data(NULL), private_Mx(my_Mx), private_My(my_My)
  { }

  ~Array2D() {
    if (wrapper == false)
      free(private_data);
  }

  void set_size(int my_Mx, int my_My) {
    private_Mx = my_Mx;
    private_My = my_My;
  }

  int allocate() {
    if (private_data != NULL)
      free(private_data);

    private_data = (T*)malloc(Mx() * My() * sizeof(T));
    if (private_data == NULL)
      return 1;

    return 0;
  }

  void wrap(T* my_data) {
    if (private_data != NULL && wrapper == false)
      free(private_data);

    private_data = my_data;
    wrapper = true;
  }

  inline T& operator()(int i, int j)
  {
    return private_data[j * Mx() + i];
  }

  inline int Mx()
  {
    return private_Mx;
  }

  inline int My()
  {
    return private_My;
  }

  inline T* data()
  {
    return private_data;
  }

private:
  bool wrapper;
  T* private_data;
  int private_Mx, private_My;
};

#endif /* _ARRAY2D_H_ */
