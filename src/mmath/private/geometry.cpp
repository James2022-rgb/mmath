// TU header --------------------------------------------
#include "mmath/public/geometry.h"

// c++ headers ------------------------------------------

// project headers --------------------------------------
#include "mmath/public/vector.h"

namespace mmath {

void MakeLookAtMatrix(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL eye,
  float const* MBASE_NOT_NULL target,
  float const* MBASE_NOT_NULL up
) {
  Vector3 const f = (Vector3(target) - Vector3(eye)).Normalize();

  MakeLookAtMatrixFromCameraOrientation(
    dst,
    eye,
    f.Data(),
    up
  );
}

void MakeLookAtMatrixFromCameraOrientation(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL eye,
  float const* MBASE_NOT_NULL forward,
  float const* MBASE_NOT_NULL up
) {
  Vector3 const e(eye);
  Vector3 const f = -Vector3(forward);
  Vector3 const u(up);
  Vector3 const l = Vector3::Cross(u, f);

  dst[0] = l.X();
  dst[1] = u.X();
  dst[2] = f.X();
  dst[3] = 0.0f;
  dst[4] = l.Y();
  dst[5] = u.Y();
  dst[6] = f.Y();
  dst[7] = 0.0f;
  dst[8] = l.Z();
  dst[9] = u.Z();
  dst[10] = f.Z();
  dst[11] = 0.0f;
  dst[12] = -Vector3::Dot(l, e);
  dst[13] = -Vector3::Dot(u, e);
  dst[14] = -Vector3::Dot(f, e);
  dst[15] = 1.0f;
}

} // namespace mmath
