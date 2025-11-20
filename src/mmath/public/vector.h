#pragma once

// c++ headers ------------------------------------------
#include <array>

// public project headers -------------------------------
#include "mbase/public/type_util.h"

namespace mmath {

struct Vector2 {
  // Technically UB, but widely supported.
  union {
    struct { float x, y; };
    std::array<float, 2> data;
  };

  Vector2() : data{0.0f, 0.0f} {}
  Vector2(float x_, float y_) : data{x_, y_} {}
  Vector2(float const* MBASE_NOT_NULL v) : data{v[0], v[1]} {}

  Vector2(Vector2 const& other) : data(other.data) {}
  Vector2& operator=(Vector2 const& other) {
    if (this != &other) {
      data = other.data;
    }
    return *this;
  }

  Vector2 operator-() const {
    return Vector2(-data[0], -data[1]);
  }

  Vector2 operator+(Vector2 const& other) const {
    return Vector2(data[0] + other.data[0], data[1] + other.data[1]);
  }
  Vector2 operator-(Vector2 const& other) const {
    return Vector2(data[0] - other.data[0], data[1] - other.data[1]);
  }
  Vector2 operator*(Vector2 const& other) const {
    return Vector2(data[0] * other.data[0], data[1] * other.data[1]);
  }
  Vector2 operator/(Vector2 const& other) const {
    return Vector2(data[0] / other.data[0], data[1] / other.data[1]);
  }

  Vector2 operator*(float scalar) const {
    return Vector2(data[0] * scalar, data[1] * scalar);
  }
  Vector2 operator/(float scalar) const {
    return Vector2(data[0] / scalar, data[1] / scalar);
  }

  Vector2& operator+=(Vector2 const& other) {
    data[0] += other.data[0];
    data[1] += other.data[1];
    return *this;
  }
  Vector2& operator-=(Vector2 const& other) {
    data[0] -= other.data[0];
    data[1] -= other.data[1];
    return *this;
  }

  Vector2& operator*=(float scalar) {
    data[0] *= scalar;
    data[1] *= scalar;
    return *this;
  }
  Vector2& operator/=(float scalar) {
    data[0] /= scalar;
    data[1] /= scalar;
    return *this;
  }

  inline float& operator[](size_t index) {
    return data[index];
  }

  inline float const& operator[](size_t index) const {
    return data[index];
  }

  inline float& X() { return data[0]; }
  inline float& Y() { return data[1]; }

  inline float const& X() const { return data[0]; }
  inline float const& Y() const { return data[1]; }

  inline float const* MBASE_NOT_NULL Data() const {
    return data.data();
  }

  float Length2() const {
    return data[0] * data[0] + data[1] * data[1];
  }
  float Length() const;

  Vector2 Normalize() const {
    float len = this->Length();
    if (len > 0.0f) {
      return Vector2(data[0] / len, data[1] / len);
    }
    return *this;
  }

  static float Dot(
    Vector2 const& a,
    Vector2 const& b
  ) {
    return
      a.data[0] * b.data[0] +
      a.data[1] * b.data[1];
  }

  static float Cross(
    Vector2 const& a,
    Vector2 const& b
  ) {
    return a.data[0] * b.data[1] - a.data[1] * b.data[0];
  }
};

static_assert(sizeof(Vector2) == 2 * sizeof(float));
static_assert(alignof(Vector2) == alignof(std::array<float, 2>));

struct Vector3 {
  // Technically UB, but widely supported.
  union {
    struct { float x, y, z; };
    std::array<float, 3> data;
  };

  Vector3() : data{0.0f, 0.0f, 0.0f} {}
  Vector3(float x_, float y_, float z_) : data{x_, y_, z_} {}
  Vector3(float const* MBASE_NOT_NULL v) : data{v[0], v[1], v[2]} {}

  Vector3(Vector3 const& other) : data(other.data) {}
  Vector3& operator=(Vector3 const& other) {
    if (this != &other) {
      data = other.data;
    }
    return *this;
  }

  Vector3 operator-() const {
    return Vector3(-data[0], -data[1], -data[2]);
  }

  Vector3 operator+(Vector3 const& other) const {
    return Vector3(data[0] + other.data[0], data[1] + other.data[1], data[2] + other.data[2]);
  }
  Vector3 operator-(Vector3 const& other) const {
    return Vector3(data[0] - other.data[0], data[1] - other.data[1], data[2] - other.data[2]);
  }
  Vector3 operator*(Vector3 const& other) const {
    return Vector3(data[0] * other.data[0], data[1] * other.data[1], data[2] * other.data[2]);
  }
  Vector3 operator/(Vector3 const& other) const {
    return Vector3(data[0] / other.data[0], data[1] / other.data[1], data[2] / other.data[2]);
  }

  Vector3 operator*(float scalar) const {
    return Vector3(data[0] * scalar, data[1] * scalar, data[2] * scalar);
  }
  Vector3 operator/(float scalar) const {
    return Vector3(data[0] / scalar, data[1] / scalar, data[2] / scalar);
  }

  Vector3& operator+=(Vector3 const& other) {
    data[0] += other.data[0];
    data[1] += other.data[1];
    data[2] += other.data[2];
    return *this;
  }
  Vector3& operator-=(Vector3 const& other) {
    data[0] -= other.data[0];
    data[1] -= other.data[1];
    data[2] -= other.data[2];
    return *this;
  }

  Vector3& operator*=(float scalar) {
    data[0] *= scalar;
    data[1] *= scalar;
    data[2] *= scalar;
    return *this;
  }
  Vector3& operator/=(float scalar) {
    data[0] /= scalar;
    data[1] /= scalar;
    data[2] /= scalar;
    return *this;
  }

  inline float& operator[](size_t index) {
    return data[index];
  }

  inline float const& operator[](size_t index) const {
    return data[index];
  }

  inline float& X() { return data[0]; }
  inline float& Y() { return data[1]; }
  inline float& Z() { return data[2]; }

  inline float const& X() const { return data[0]; }
  inline float const& Y() const { return data[1]; }
  inline float const& Z() const { return data[2]; }

  inline float* MBASE_NOT_NULL Data() {
    return data.data();
  }
  inline float const* MBASE_NOT_NULL Data() const {
    return data.data();
  }

  float Length2() const {
    return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
  }
  float Length() const;

  Vector3 Normalize() const {
    float len = this->Length();
    if (len > 0.0f) {
      return Vector3(data[0] / len, data[1] / len, data[2] / len);
    }
    return *this;
  }

  static float Dot(
    Vector3 const& a,
    Vector3 const& b
  ) {
    return
      a.data[0] * b.data[0] +
      a.data[1] * b.data[1] +
      a.data[2] * b.data[2];
  }

  static Vector3 Cross(
    Vector3 const& a,
    Vector3 const& b
  ) {
    return Vector3(
      a.data[1] * b.data[2] - a.data[2] * b.data[1],
      a.data[2] * b.data[0] - a.data[0] * b.data[2],
      a.data[0] * b.data[1] - a.data[1] * b.data[0]
    );
  }
};

static_assert(sizeof(Vector3) == 3 * sizeof(float));
static_assert(alignof(Vector3) == alignof(std::array<float, 3>));

struct Vector4 {
  // Technically UB, but widely supported.
  union {
    struct { float x, y, z, w; };
    std::array<float, 4> data;
  };

  Vector4() : data{0.0f, 0.0f, 0.0f, 0.0f} {}
  Vector4(float x_, float y_, float z_, float w_) : data{x_, y_, z_, w_} {}
  Vector4(float const* MBASE_NOT_NULL v) : data{v[0], v[1], v[2], v[3]} {}
  Vector4(Vector3 const& v3, float w_) : data{ v3.data[0], v3.data[1], v3.data[2], w_ } {}

  Vector4(Vector4 const& other) : data(other.data) {}
  Vector4& operator=(Vector4 const& other) {
    if (this != &other) {
      data = other.data;
    }
    return *this;
  }

  Vector4 operator-() const {
    return Vector4(-data[0], -data[1], -data[2], -data[3]);
  }

  Vector4 operator+(Vector4 const& other) const {
    return Vector4(data[0] + other.data[0], data[1] + other.data[1], data[2] + other.data[2], data[3] + other.data[3]);
  }
  Vector4 operator-(Vector4 const& other) const {
    return Vector4(data[0] - other.data[0], data[1] - other.data[1], data[2] - other.data[2], data[3] - other.data[3]);
  }
  Vector4 operator*(Vector4 const& other) const {
    return Vector4(data[0] * other.data[0], data[1] * other.data[1], data[2] * other.data[2], data[3] * other.data[3]);
  }
  Vector4 operator/(Vector4 const& other) const {
    return Vector4(data[0] / other.data[0], data[1] / other.data[1], data[2] / other.data[2], data[3] / other.data[3]);
  }

  Vector4 operator*(float scalar) const {
    return Vector4(data[0] * scalar, data[1] * scalar, data[2] * scalar, data[3] * scalar);
  }
  Vector4 operator/(float scalar) const {
    return Vector4(data[0] / scalar, data[1] / scalar, data[2] / scalar, data[3] / scalar);
  }

  Vector4& operator+=(Vector4 const& other) {
    data[0] += other.data[0];
    data[1] += other.data[1];
    data[2] += other.data[2];
    data[3] += other.data[3];
    return *this;
  }
  Vector4& operator-=(Vector4 const& other) {
    data[0] -= other.data[0];
    data[1] -= other.data[1];
    data[2] -= other.data[2];
    data[3] -= other.data[3];
    return *this;
  }

  Vector4& operator*=(float scalar) {
    data[0] *= scalar;
    data[1] *= scalar;
    data[2] *= scalar;
    data[3] *= scalar;
    return *this;
  }
  Vector4& operator/=(float scalar) {
    data[0] /= scalar;
    data[1] /= scalar;
    data[2] /= scalar;
    data[3] /= scalar;
    return *this;
  }

  inline float& operator[](size_t index) {
    return data[index];
  }

  inline float const& operator[](size_t index) const {
    return data[index];
  }

  inline float& X() { return data[0]; }
  inline float& Y() { return data[1]; }
  inline float& Z() { return data[2]; }

  inline float const& X() const { return data[0]; }
  inline float const& Y() const { return data[1]; }
  inline float const& Z() const { return data[2]; }

  inline float* MBASE_NOT_NULL Data() {
    return data.data();
  }
  inline float const* MBASE_NOT_NULL Data() const {
    return data.data();
  }

  float Length2() const {
    return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
  }
  float Length() const;

  Vector4 Normalize() const {
    float len = this->Length();
    if (len > 0.0f) {
      return Vector4(data[0] / len, data[1] / len, data[2] / len, data[3] / len);
    }
    return *this;
  }

  static float Dot(
    Vector4 const& a,
    Vector4 const& b
  ) {
    return
      a.data[0] * b.data[0] +
      a.data[1] * b.data[1] +
      a.data[2] * b.data[2] +
      a.data[3] * b.data[3];
  }
};

static_assert(sizeof(Vector4) == 4 * sizeof(float));
static_assert(alignof(Vector4) == alignof(std::array<float, 4>));

} // namespace mmath
