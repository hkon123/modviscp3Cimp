
import numpy as np
cimport numpy as np

cdef extern from "test.c":
   float c_f "f"(float *, float *, float *, int)
   int c_j "j"(float *, float *, float *,float *, float *, float)
   int c_g "g"(float *, float *, float *,float *, float *, float)

cpdef f(in_phi, iterations):
  cdef np.ndarray[np.float32_t,ndim=2] out_phi = np.empty((100, 100), dtype=np.float32)
  cdef np.ndarray[np.float32_t,ndim=2] np_in_phi = np.ascontiguousarray(in_phi, dtype=np.float32)
  cdef np.ndarray[np.float32_t,ndim=2] out_en = np.empty((100, 100), dtype=np.float32)

  cdef float energysum = c_f(<float *> np_in_phi.data, <float *> out_phi.data, <float *> out_en.data, iterations)

  return out_phi, out_en, energysum

cpdef j(in_rhoj, treshold):
  cdef np.ndarray[np.float32_t,ndim=3] out_phij = np.empty((50, 50, 50), dtype=np.float32)
  cdef np.ndarray[np.float32_t,ndim=1] out_conv = np.empty((10000), dtype=np.float32)
  cdef np.ndarray[np.float32_t,ndim=3] np_in_rhoj = np.ascontiguousarray(in_rhoj, dtype=np.float32)
  cdef np.ndarray[np.float32_t,ndim=3] out_Exj = np.empty((50, 50, 50), dtype=np.float32)
  cdef np.ndarray[np.float32_t,ndim=3] out_Eyj = np.empty((50, 50, 50), dtype=np.float32)

  cdef int value = c_j(<float *> np_in_rhoj.data,<float *> out_phij.data, <float *> out_conv.data, <float *> out_Exj.data, <float *> out_Eyj.data, treshold)

  return out_phij, out_conv, out_Exj, out_Eyj, value

cpdef g(in_rhog, tresholdg):
  cdef np.ndarray[np.float32_t,ndim=3] out_phig = np.empty((50, 50, 50), dtype=np.float32)
  cdef np.ndarray[np.float32_t,ndim=1] out_convg = np.empty((10000), dtype=np.float32)
  cdef np.ndarray[np.float32_t,ndim=3] np_in_rhog = np.ascontiguousarray(in_rhog, dtype=np.float32)
  cdef np.ndarray[np.float32_t,ndim=3] out_Exg = np.empty((50, 50, 50), dtype=np.float32)
  cdef np.ndarray[np.float32_t,ndim=3] out_Eyg = np.empty((50, 50, 50), dtype=np.float32)

  cdef int value = c_g(<float *> np_in_rhog.data, <float *> out_phig.data, <float *> out_convg.data, <float *> out_Exg.data, <float *> out_Eyg.data, tresholdg)

  return out_phig, out_convg, out_Exg, out_Eyg, value
