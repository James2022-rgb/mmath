// TU header --------------------------------------------
#include "mmath/public/projection.h"

// c++ headers ------------------------------------------
#include <cmath>

// compile-time flags -----------------------------------

//
// The Perspective and Orthographic Projection Matrix
// https://www.scratchapixel.com/lessons/3d-basic-rendering/perspective-and-orthographic-projection-matrix/building-basic-perspective-projection-matrix.html
// Projection Matrixについて | shikihuiku – 色不異空 – Real-time rendering topics in Japanese.
// https://shikihuiku.github.io/post/projection_matrix/
//

namespace mmath {

void MakePerspectiveProjectionMatrixZeroOneRH(float* MBASE_NOT_NULL dst, float aspect_ratio, float fov_v, float near_plane, float far_plane) {
  //
  // For depth range [-1, 1]:
  //
  // | 1 / (aspect_ratio * tan(fov / 2)),          0,                     0,                               0              |
  // |               0,                   1 / (tan(fov / 2),              0,                               0              |
  // |               0,                            0,        (-near -far) / (near - far), (2 * far * near) / (near - far) |
  // |               0,                            0,                     1,                               0              |
  //
  // For depth range [0, 1]:
  //
  // | 1 / (aspect_ratio * tan(fov / 2)),          0,                   0,                    0              |
  // |               0,                   1 / (tan(fov / 2),            0,                    0              |
  // |               0,                            0,        far / (near - far), (far * near) / (near - far) |
  // |               0,                            0,                   1,                    0              |
  //
  // https://stackoverflow.com/questions/11277501/how-to-recover-view-space-position-given-view-space-depth-value-and-ndc-xy
  //

  float const tan_half_fov_v = std::tan(fov_v / 2.0f);
  float const near_minus_far = near_plane - far_plane;

  float const handedness = -1.0f; // Right-handed coordinate system.

  dst[0] = 1.0f / (aspect_ratio * tan_half_fov_v); // Scale the x coordinates of the projected point.
  dst[5] = 1.0f / tan_half_fov_v * -1.0f;          // Scale the y coordinates of the projected point.
  dst[10] = far_plane / near_minus_far;
  dst[11] = handedness; // Right-handed: Flip the z-axis.
  dst[14] = (far_plane * near_plane) / near_minus_far;
  dst[1] = 0.0f;
  dst[2] = 0.0f;
  dst[3] = 0.0f;
  dst[4] = 0.0f;
  dst[6] = 0.0f;
  dst[7] = 0.0f;
  dst[8] = 0.0f;
  dst[9] = 0.0f;
  dst[12] = 0.0f;
  dst[13] = 0.0f;
  dst[15] = 0.0f;
}

} // namespace mmath
