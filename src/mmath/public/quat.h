#pragma once

// c++ headers ------------------------------------------
#include <array>

// public project headers -------------------------------
#include "mbase/public/type_util.h"

// project headers --------------------------------------
#include "mmath/public/vector.h"
#include "mmath/public/matrix.h"

namespace mmath {

//
// ## Remarks
// A quaternion is represented as an array of 4 floats: [x, y, z, w] where x, y, z comprise the vector (imaginary) part and w the scalar (real) part.
//

struct Quat final {
  std::array<float, 4> data;

  Quat() : data{ 0.0f, 0.0f, 0.0f, 1.0f } {} // Default to identity quaternion
  Quat(float x, float y, float z, float w) : data{ x, y, z, w } {}
  Quat(float const* MBASE_NOT_NULL v) : data{ v[0], v[1], v[2], v[3] } {}

  Quat(Quat const& other) : data(other.data) {}
  Quat& operator=(Quat const& other) {
    if (this != &other) {
      data = other.data;
    }
    return *this;
  }

  Quat(Quat&& other) noexcept : data(std::move(other.data)) {}
  Quat& operator=(Quat&& other) noexcept {
    if (this != &other) {
      data = std::move(other.data);
    }
    return *this;
  }

  static Quat AxisAngle(
    float const* MBASE_NOT_NULL axis,
    float angle
  );
  static Quat AxisAngle(
    Vector3 const& axis,
    float angle
  ) {
    return AxisAngle(axis.Data(), angle);
  }

  Quat operator-() const {
    return Quat(-data[0], -data[1], -data[2], data[3]);
  }

  Quat operator-(Quat const& other) const {
    return Quat(data[0] - other.data[0], data[1] - other.data[1], data[2] - other.data[2], data[3] - other.data[3]);
  }
  Quat operator+(Quat const& other) const {
    return Quat(data[0] + other.data[0], data[1] + other.data[1], data[2] + other.data[2], data[3] + other.data[3]);
  }

  Quat operator*(float scalar) const {
    return Quat(data[0] * scalar, data[1] * scalar, data[2] * scalar, data[3] * scalar);
  }
  Quat operator/(float scalar) const {
    return Quat(data[0] / scalar, data[1] / scalar, data[2] / scalar, data[3] / scalar);
  }

  Quat& operator+=(Quat const& other) {
    data[0] += other.data[0];
    data[1] += other.data[1];
    data[2] += other.data[2];
    data[3] += other.data[3];
    return *this;
  }
  Quat& operator-=(Quat const& other) {
    data[0] -= other.data[0];
    data[1] -= other.data[1];
    data[2] -= other.data[2];
    data[3] -= other.data[3];
    return *this;
  }

  Quat& operator*=(float scalar) {
    data[0] *= scalar;
    data[1] *= scalar;
    data[2] *= scalar;
    data[3] *= scalar;
    return *this;
  }
  Quat& operator/=(float scalar) {
    data[0] /= scalar;
    data[1] /= scalar;
    data[2] /= scalar;
    data[3] /= scalar;
    return *this;
  }

  Quat operator*(Quat const& other) const;
  Quat& operator*=(Quat const& other);

  inline float& X() { return data[0]; }
  inline float& Y() { return data[1]; }
  inline float& Z() { return data[2]; }
  inline float& W() { return data[3]; }

  inline float const& X() const { return data[0]; }
  inline float const& Y() const { return data[1]; }
  inline float const& Z() const { return data[2]; }
  inline float const& W() const { return data[3]; }

  inline float const* MBASE_NOT_NULL Data() const {
    return data.data();
  }
  inline float* MBASE_NOT_NULL Data() {
    return data.data();
  }

  void ToMat33(float* MBASE_NOT_NULL dst) const;
  Matrix33 ToMat33() const;

  void ToMat44(float* MBASE_NOT_NULL dst) const;
  Matrix44 ToMat44() const;

  float Norm() const;

  void NormalizeInplace();
  bool Invert();

  void Normalized(Quat& dst) const;
  bool Inverted(Quat& dst) const;
};

Vector3 operator*(Quat const& q, Vector3 const& v);

void IdentityQuat(float* MBASE_NOT_NULL dst);

float NormQuat(
  float const* MBASE_NOT_NULL q
);

void NormalizeQuat(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL q
);

void ConjugateQuat(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL q
);

void MulQuatQuat(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL a,
  float const* MBASE_NOT_NULL b
);

bool InvertQuat(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL q
);

bool InvertQuatInplace(
  float* MBASE_NOT_NULL q
);

/// ## Parameters
/// - `dst`: Pointer to the destination quaternion (x, y, z, w).
/// - `axis`: Pointer to the axis of rotation (3 floats) which MUST be normalized.
/// - `angle`: The angle of rotation in radians.
void RotateAxisAngleQuat(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL axis,
  float angle
);

void QuatToMat33(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL q
);
void QuatToMat44(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL q
);

void ApplyQuatToVec3(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL q,
  float const* MBASE_NOT_NULL v
);

} // namespace mmath

