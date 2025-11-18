#pragma once

// public project headers -------------------------------
#include "mbase/public/type_util.h"

// project headers --------------------------------------
#include "mmath/public/vector.h"
#include "mmath/public/quat.h"
#include "mmath/public/matrix.h"

namespace mmath {

inline void RotateAxisAngleQuat(
  float* MBASE_NOT_NULL dst,
  Vector3 const& axis,
  float angle
) {
  RotateAxisAngleQuat(
    dst,
    axis.Data(),
    angle
  );
}

/// Creates a look-at matrix.
/// The resulting matrix transforms a point in the world space to the camera space,
/// such that the camera is positioned at `eye`, and its -Z axis points towards `target`.
///
/// ## Parameters
/// - `dst`: Pointer to the destination matrix (16 floats).
/// - `eye`: Pointer to the eye position (3 floats).
/// - `target`: Pointer to the target position (3 floats).
/// - `up`: Pointer to the up vector (3 floats) which MUST be normalized.
///
/// ## Remarks
/// - The resulting matrix is multiplied to the left of a column vector, and is laid out in column-major order.
void MakeLookAtMatrix(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL eye,
  float const* MBASE_NOT_NULL target,
  float const* MBASE_NOT_NULL up
);

/// Creates a look-at matrix.
/// The resulting matrix transforms a point in the world space to the camera space,
/// such that the camera is positioned at `eye`, and its -Z axis is aligned with `forward`.
///
/// ## Parameters
/// - `dst`: Pointer to the destination matrix (16 floats).
/// - `eye`: Pointer to the eye position (3 floats).
/// - `forward`: Pointer to the forward vector (3 floats) which MUST be normalized.
/// - `up`: Pointer to the up vector (3 floats) which MUST be normalized.
///
/// ## Remarks
/// - The resulting matrix is multiplied to the left of a column vector, and is laid out in column-major order.
void MakeLookAtMatrixFromCameraOrientation(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL eye,
  float const* MBASE_NOT_NULL forward,
  float const* MBASE_NOT_NULL up
);

inline void MakeLookAtMatrix(
  float* MBASE_NOT_NULL dst,
  Vector3 const& eye,
  Vector3 const& target,
  Vector3 const& up
) {
  MakeLookAtMatrix(dst, eye.Data(), target.Data(), up.Data());
}

inline void MakeLookAtMatrixFromCameraOrientation(
  float* MBASE_NOT_NULL dst,
  Vector3 const& eye,
  Vector3 const& forward,
  Vector3 const& up
) {
  MakeLookAtMatrixFromCameraOrientation(dst, eye.Data(), forward.Data(), up.Data());
}

} // namespace mmath
