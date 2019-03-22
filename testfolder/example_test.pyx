
import numpy as np
cimport numpy as np

cdef extern from "test.c":
   float c_f "f"(float *, float *, int)

cpdef f(in_phi, iterations):
  cdef np.ndarray[np.float32_t,ndim=2] out_phi = np.empty((100, 100), dtype=np.float32)
  cdef np.ndarray[np.float32_t,ndim=2] np_in_phi = np.ascontiguousarray(in_phi, dtype=np.float32)

  cdef float energysum = c_f(<float *> np_in_phi.data, <float *> out_phi.data, iterations)

  return out_phi, energysum
