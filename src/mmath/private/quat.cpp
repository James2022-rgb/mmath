// TU header --------------------------------------------
#include "mmath/public/quat.h"

// c++ headers ------------------------------------------
#include <cstring>
#include <cmath>

//
// (Yet another) Introduction to quaternions | lisyarus blog
// https://lisyarus.github.io/blog/posts/introduction-to-quaternions.html
// クォータニオン (Quaternion) を総整理！ ～ 三次元物体の回転と姿勢を鮮やかに扱う ～ #Unity - Qiita
// https://qiita.com/drken/items/0639cf34cce14e8d58a5
//

namespace mmath {

Quat Quat::AxisAngle(
  float const* MBASE_NOT_NULL axis,
  float angle
) {
  Quat result;
  RotateAxisAngleQuat(result.data.data(), axis, angle);
  return result;
}

Quat Quat::operator*(Quat const& other) const {
  std::array<float, 4> result;
  MulQuatQuat(result.data(), this->data.data(), other.data.data());
  return Quat(result.data());
}

Quat& Quat::operator*=(Quat const& other) {
  std::array<float, 4> result;
  MulQuatQuat(result.data(), this->data.data(), other.data.data());
  this->data = result;
  return *this;
}

void Quat::ToMat33(float* MBASE_NOT_NULL dst) const {
  QuatToMat33(dst, this->data.data());
}

Matrix33 Quat::ToMat33() const {
  Matrix33 result;
  this->ToMat33(result.data.data());
  return result;
}

void Quat::ToMat44(float* MBASE_NOT_NULL dst) const {
  QuatToMat44(dst, this->data.data());
}

Matrix44 Quat::ToMat44() const {
  Matrix44 result;
  this->ToMat44(result.data.data());
  return result;
}

float Quat::Norm() const {
  return NormQuat(this->data.data());
}

void Quat::NormalizeInplace() {
  std::array<float, 4> normalized;
  NormalizeQuat(normalized.data(), this->data.data());
  this->data = normalized;
}

bool Quat::Invert() {
  std::array<float, 4> inverted;
  if (!InvertQuat(inverted.data(), this->data.data())) {
    return false; // Inversion failed (e.g., norm is zero)
  }
  this->data = inverted;
  return true;
}

void Quat::Normalized(Quat& dst) const {
  dst = *this;
  dst.NormalizeInplace();
}

bool Quat::Inverted(Quat& dst) const {
  Quat result = *this;
  if (!InvertQuatInplace(result.data.data())) {
    return false;
  }
  dst = result;
  return true;
}

Vector3 operator*(Quat const& q, Vector3 const& v) {
  Vector3 result;
  ApplyQuatToVec3(result.Data(), q.data.data(), v.Data());
  return result;
}

void IdentityQuat(float* MBASE_NOT_NULL dst) {
  dst[0] = 0.0f; // x
  dst[1] = 0.0f; // y
  dst[2] = 0.0f; // z
  dst[3] = 1.0f; // w
}

float NormQuat(
  float const* MBASE_NOT_NULL q
) {
  // Norm = sqrt(x^2 + y^2 + z^2 + w^2)
  return std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
}

void NormalizeQuat(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL q
) {
  // Normalize: q' = q / ||q||
  float const norm = NormQuat(q);
  if (norm == 0.0f) {
    IdentityQuat(dst); // Avoid division by zero
    return;
  }
  float const inv_norm = 1.0f / norm;
  dst[0] = q[0] * inv_norm; // x
  dst[1] = q[1] * inv_norm; // y
  dst[2] = q[2] * inv_norm; // z
  dst[3] = q[3] * inv_norm; // w
}

void ConjugateQuat(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL q
) {
  // Conjugate: [x, y, z, w] -> [-x, -y, -z, w]
  dst[0] = -q[0]; // x
  dst[1] = -q[1]; // y
  dst[2] = -q[2]; // z
  dst[3] =  q[3]; // w
}

void MulQuatQuat(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL a,
  float const* MBASE_NOT_NULL b
) {
  // ij =  k, jk =  i, ki =  j
  // ji = -k, kj = -i, ik = -j
  // ii = jj = kk = -1
  // ijk = -1
  //
  //   [ai + bj + ck + d] * [ei + fj + gk + h]
  // = (ah + bg - cf + de)i + (-ag + bh + ce + df)j + (af - be + ch + dg)k + (-ae - bf - cg + dh)

  dst[0] =  a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]; // x
  dst[1] = -a[0] * b[2] + a[1] * b[3] + a[2] * b[0] + a[3] * b[1]; // y
  dst[2] =  a[0] * b[1] - a[1] * b[0] + a[2] * b[3] + a[3] * b[2]; // z
  dst[3] = -a[0] * b[0] - a[1] * b[1] - a[2] * b[2] + a[3] * b[3]; // w
}

bool InvertQuat(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL q
) {
  memcpy(dst, q, 4 * sizeof(float)); // Copy the quaternion to the destination
  return InvertQuatInplace(dst); // Invert in place
}

bool InvertQuatInplace(
  float* MBASE_NOT_NULL q
) {
  // Inverse: q^-1 = q* / ||q||^2
  // where q* is the conjugate of q, and ||q||^2 is the norm squared.
  float const norm = NormQuat(q);
  if (std::abs(norm) < 1e-6f) {
    return false;
  }
  float const inv_norm_squared = 1.0f / (norm * norm);
  
  q[0] = -q[0] * inv_norm_squared; // x
  q[1] = -q[1] * inv_norm_squared; // y
  q[2] = -q[2] * inv_norm_squared; // z
  q[3] =  q[3] * inv_norm_squared; // w

  return true;
}

void RotateAxisAngleQuat(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL axis,
  float angle
) {
  float const half_angle = angle * 0.5f;
  float const s = std::sin(half_angle);
  float const c = std::cos(half_angle);

  dst[0] = axis[0] * s; // x
  dst[1] = axis[1] * s; // y
  dst[2] = axis[2] * s; // z
  dst[3] = c;           // w
}

void QuatToMat33(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL q
) {
  float const x = q[0];
  float const y = q[1];
  float const z = q[2];
  float const w = q[3];

  dst[0] = 1.0f - 2.0f * (y * y + z * z);
  dst[1] = 2.0f * (x * y - w * z);
  dst[2] = 2.0f * (x * z + w * y);
  
  dst[3] = 2.0f * (x * y + w * z);
  dst[4] = 1.0f - 2.0f * (x * x + z * z);
  dst[5] = 2.0f * (y * z - w * x);

  dst[6] = 2.0f * (x * z - w * y);
  dst[7] = 2.0f * (y * z + w * x);
  dst[8] = 1.0f - 2.0f * (x * x + y * y);
}

void QuatToMat44(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL q
) {
  float const x = q[0];
  float const y = q[1];
  float const z = q[2];
  float const w = q[3];

  dst[0] = 1.0f - 2.0f * (y * y + z * z);
  dst[1] = 2.0f * (x * y - w * z);
  dst[2] = 2.0f * (x * z + w * y);
  
  dst[4] = 2.0f * (x * y + w * z);
  dst[5] = 1.0f - 2.0f * (x * x + z * z);
  dst[6] = 2.0f * (y * z - w * x);

  dst[8] = 2.0f * (x * z - w * y);
  dst[9] = 2.0f * (y * z + w * x);
  dst[10] = 1.0f - 2.0f * (x * x + y * y);

  dst[3] = 0.0f;
  dst[7] = 0.0f;

  dst[11] = 0.0f;
  dst[12] = 0.0f;
  dst[13] = 0.0f;
  dst[14] = 0.0f;
  dst[15] = 1.0f;
}

void ApplyQuatToVec3(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL q,
  float const* MBASE_NOT_NULL v
) {
  // https://gamedev.stackexchange.com/questions/28395/rotating-vector3-by-a-quaternion

  std::array<float, 4> p;
  p[0] = v[0]; // x
  p[1] = v[1]; // y
  p[2] = v[2]; // z
  p[3] = 0.0f; // w (scalar part is zero for vector)

  // q* (conjugate of q)
  std::array<float, 4> q_conjugate;
  ConjugateQuat(q_conjugate.data(), q); 

  // p' = q x p x q*

  std::array<float, 4> tmp, tmp2;
  MulQuatQuat(tmp.data(), q_conjugate.data(), p.data());
  MulQuatQuat(tmp2.data(), tmp.data(), q);

  dst[0] = tmp2[0]; // x
  dst[1] = tmp2[1]; // y
  dst[2] = tmp2[2]; // z
}

} // namespace mmath
