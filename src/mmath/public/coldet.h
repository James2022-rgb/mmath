#pragma once

// c++ headers ------------------------------------------
#include <array>

// public project headers -------------------------------
#include "mbase/public/type_util.h"

#include "mmath/public/vector.h"

namespace mmath {

struct IntersectRaySphereResult {
  Vector3 point;
  Vector3 normal;
  float t;
};

/// ## Parameters
/// - `ray_origin`: Origin of the ray (3 floats).
/// - `ray_direction`: Normalized direction of the ray (3 floats).
/// - `sphere_center`: Center of the sphere (3 floats).
///
/// ## Returns
/// - `true` if the ray intersects the sphere, `false` otherwise.
bool IntersectRaySphere(
  float const* MBASE_NOT_NULL ray_origin,
  float const* MBASE_NOT_NULL ray_direction,
  float const* MBASE_NOT_NULL sphere_center,
  float sphere_radius,
  IntersectRaySphereResult* MBASE_NULLABLE out_result
);

bool IntersectRaySphere(
  Vector3 const& ray_origin,
  Vector3 const& ray_direction,
  Vector3 const& sphere_center,
  float sphere_radius,
  IntersectRaySphereResult* MBASE_NULLABLE out_result
);

} // namespace mmath
