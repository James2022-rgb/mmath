#pragma once

// public project headers -------------------------------
#include "mbase/public/type_util.h"

namespace mmath {

/// ## Remarks
/// - Right-handed NDC space is used.
/// - Z values will be in the range [0, 1] after perspective division.
/// - The resulting matrix is multiplied to the left of a column vector, and is laid out in column-major order.
void MakePerspectiveProjectionMatrixZeroOneRH(float* MBASE_NOT_NULL dst, float aspect_ratio, float fov_v, float near_plane, float far_plane);

} // namespace mmath
