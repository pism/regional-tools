#include <cstdlib>

#ifndef _ARRAY2D_H_
#define _ARRAY2D_H_

/* A zero-cost 2D array class hiding storage order details. */
/* This is a way to have clean code and avoid using evil macros. */
template<class T>
class Array2D {
public:
  Array2D(int Mx, int My)
    : m_wrapper(false), m_data(NULL), m_Mx(Mx), m_My(My)
  { }

  ~Array2D() {
    if (not m_wrapper) {
      free(m_data);
    }
  }

  void set_size(int Mx, int My) {
    m_Mx = Mx;
    m_My = My;
  }

  int allocate() {
    if (m_data != NULL)
      free(m_data);

    m_data = (T*)malloc(Mx() * My() * sizeof(T));
    if (m_data == NULL)
      return 1;

    return 0;
  }

  void wrap(T* data) {
    if (m_data != NULL && not m_wrapper)
      free(m_data);

    m_data = data;
    m_wrapper = true;
  }

  inline T& operator()(int i, int j) {
    return m_data[j * Mx() + i];
  }

  inline int Mx() {
    return m_Mx;
  }

  inline int My() {
    return m_My;
  }

  inline T* data() {
    return m_data;
  }

private:
  bool m_wrapper;
  T* m_data;
  int m_Mx, m_My;
};

#endif /* _ARRAY2D_H_ */
