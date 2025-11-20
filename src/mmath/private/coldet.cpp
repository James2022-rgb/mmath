// TU header --------------------------------------------
#include "mmath/public/coldet.h"

// c++ headers ------------------------------------------
#include <cmath>

// public project headers -------------------------------
#include "mmath/public/vector.h"

namespace mmath {

bool IntersectRaySphere(
  float const* MBASE_NOT_NULL ray_origin,
  float const* MBASE_NOT_NULL _ray_direction,
  float const* MBASE_NOT_NULL sphere_center,
  float sphere_radius,
  IntersectRaySphereResult* MBASE_NULLABLE out_result
) {
  Vector3 ray_direction(_ray_direction);
  if (ray_direction.Length2() <= 1e-12f) {
    // Degenerate ray direction.
    return false;
  }

  Vector3 const A(ray_origin);
  Vector3 const B = A + Vector3(ray_direction);
  Vector3 const C(sphere_center);
  float const r = sphere_radius;

  //
  // Equation of a sphere:
  //  (x-C.x)^2 + (y-C.y)^2 + (z-C.z)^2 = r^2
  // 
  // Parameteric equations for the 3D ray (where t is the parameter):
  //  x = A.x + t(B.x - A.x)
  //  y = A.y + t(B.y - A.y)
  //  z = A.z + t(B.z - A.z)
  // 
  // Substituting the ray equations into the sphere equation gives us a quadratic equation in t:
  //  (A.x + t(B.x - A.x) - C.x)^2 + (A.y + t(B.y - A.y) - C.y)^2 + (A.z + t(B.z - A.z) - C.z)^2 = r^2
  // 
  // Substituting D = B - A, E = A - C for clarity:
  //  (A.x + tD.x - C.x)^2 + (A.y + tD.y - C.y)^2 + (A.z + tD.z - C.z)^2 = r^2
  //  (E.x + tD.x)^2 + (E.y + tD.y)^2 + (E.z + tD.z)^2 = r^2
  // 
  // Expanding and rearranging gives:
  //  (D.x^2 + D.y^2 + D.z^2)t^2 + 2(D.xE.x + D.yE.y + D.zE.z)t + (E.x^2 + E.y^2 + E.z^2 - r^2) = 0
  // 
  // Discriminant (d) of the quadratic equation:
  //  d = b^2 - 4ac
  //
  // If d < 0, there is no intersection.
  // If d = 0, there is one intersection point (the ray is tangent to the sphere).
  // If d > 0, there are two intersection points.
  //

  Vector3 const D = B - A;
  Vector3 const E = A - C;

  // Coefficients of the quadratic equation:
  float const a = 1.0f; // D.x * D.x + D.y * D.y + D.z * D.z ; We assume that the direction vector is normalized.
  float const b = 2.0f * Vector3::Dot(D, E);
  float const c = Vector3::Dot(E, E) - r * r;

  // Discriminant
  float const delta = b * b - 4.0f * a * c;

  // We treat tangent as no intersection for rays.
  if (delta <= 0) {
    return false;
  }

  float const delta2 = std::sqrt(delta);

  // Compute the two intersection distances along the ray.
  float const t1 = (-b - delta2) / (2.0f * a);
  float const t2 = (-b + delta2) / (2.0f * a);

  float t;
  if (0.0f <= t1) {
    // The first intersection is in front of the ray origin.
    t = t1;
  }
  else if (0.0f <= t2) {
    // The second intersection is in front of the ray origin.
    // The ray origin is inside the sphere.
    return false;
  }
  else {
    // Both intersections are behind the ray origin.
    return false;
  }

  if (out_result != nullptr) {
    Vector3 p = A + (B - A) * t;
    Vector3 n = (p - C).Normalize();
    out_result->point = p;
    out_result->normal = n;
    out_result->t = t;
  }

  return true;
}

bool IntersectRaySphere(
  Vector3 const& ray_origin,
  Vector3 const& ray_direction,
  Vector3 const& sphere_center,
  float sphere_radius,
  IntersectRaySphereResult* MBASE_NULLABLE out_result
) {
  return IntersectRaySphere(
    ray_origin.Data(),
    ray_direction.Data(),
    sphere_center.Data(),
    sphere_radius,
    out_result
  );
}

} // namespace mmath
