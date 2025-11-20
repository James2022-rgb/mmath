#pragma once

// c++ headers ------------------------------------------
#include <array>

// public project headers -------------------------------
#include "mbase/public/type_util.h"

#include "mmath/public/vector.h"

namespace mmath {

struct Matrix33 final {
  std::array<float, 9> data;

  Matrix33() {
    data.fill(0.0f);
    data[0] = 1.0f; data[4] = 1.0f; data[8] = 1.0f;
  }
  ~Matrix33() = default;

  Matrix33(Matrix33 const& other) : data(other.data) {}
  Matrix33& operator=(Matrix33 const& other) {
    if (this != &other) {
      data = other.data;
    }
    return *this;
  }

  Matrix33(Matrix33&& other) noexcept : data(std::move(other.data)) {}
  Matrix33& operator=(Matrix33&& other) noexcept {
    if (this != &other) {
      data = std::move(other.data);
    }
    return *this;
  }

  Matrix33(float const* MBASE_NOT_NULL v);
  Matrix33(
    float const* MBASE_NOT_NULL col0,
    float const* MBASE_NOT_NULL col1,
    float const* MBASE_NOT_NULL col2
  );

  float* MBASE_NOT_NULL Data() {
    return data.data();
  }
  float const* MBASE_NOT_NULL Data() const {
    return data.data();
  }

  float& operator[](std::size_t index) {
    return data[index];
  }
  float const& operator[](std::size_t index) const {
    return data[index];
  }

  Matrix33 Transpose() const;
  bool Invert(Matrix33& out_result) const;

  Matrix33 operator*(float scalar) const;
  Matrix33 operator*(Matrix33 const& other) const;

  friend Vector3 operator*(Matrix33 const& mat, Vector3 const& vec);
};

static_assert(sizeof(Matrix33) == 9 * sizeof(float));
static_assert(alignof(Matrix33) == alignof(std::array<float, 9>));

struct Matrix44 final {
  std::array<float, 16> data;

  Matrix44() {
    data.fill(0.0f);
    data[0] = 1.0f; data[5] = 1.0f; data[10] = 1.0f; data[15] = 1.0f;;
  }
  ~Matrix44() = default;

  Matrix44(Matrix44 const& other) : data(other.data) {}
  Matrix44& operator=(Matrix44 const& other) {
    if (this != &other) {
      data = other.data;
    }
    return *this;
  }

  Matrix44(Matrix44&& other) noexcept : data(std::move(other.data)) {}
  Matrix44& operator=(Matrix44&& other) noexcept {
    if (this != &other) {
      data = std::move(other.data);
    }
    return *this;
  }

  Matrix44(float const* MBASE_NOT_NULL v);
  Matrix44(
    float const* MBASE_NOT_NULL col0,
    float const* MBASE_NOT_NULL col1,
    float const* MBASE_NOT_NULL col2,
    float const* MBASE_NOT_NULL col3
  );

  float* MBASE_NOT_NULL Data() {
    return data.data();
  }
  float const* MBASE_NOT_NULL Data() const {
    return data.data();
  }

  float& operator[](std::size_t index) {
    return data[index];
  }
  float const& operator[](std::size_t index) const {
    return data[index];
  }

  Matrix44 Transpose() const;
  bool Invert(Matrix44& out_result) const;

  Matrix44 operator*(float scalar) const;
  Matrix44 operator*(Matrix44 const& other) const;

  friend Vector4 operator*(Matrix44 const& mat, Vector4 const& vec);
};

static_assert(sizeof(Matrix44) == 16 * sizeof(float));
static_assert(alignof(Matrix44) == alignof(std::array<float, 16>));

void IdentityMat33(float* MBASE_NOT_NULL dst);
void IdentityMat44(float* MBASE_NOT_NULL dst);

/// Transposes a 4x4 matrix.
///
/// ## Parameters
/// - `dst`: Pointer to the destination matrix (16 floats).
void TransposeMat44(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL src
);

/// Transposes a 3x3 matrix in place.
///
/// ## Parameters
/// - `dst`: Pointer to the matrix to transpose (9 floats).
void TransposeMat33Inplace(float* MBASE_NOT_NULL dst);

/// Transposes a 4x4 matrix in place.
///
/// ## Parameters
/// - `dst`: Pointer to the matrix to transpose (16 floats).
void TransposeMat44Inplace(float* MBASE_NOT_NULL dst);

void MulMat33Vec3(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL mat,
  float const* MBASE_NOT_NULL vec
);

void MulMat44Vec4(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL mat,
  float const* MBASE_NOT_NULL vec
);

/// ## Parameters
/// - `dst`: Pointer to the destination matrix (9 floats).
/// - `a`: Pointer to the first source matrix (9 floats).
/// - `b`: Pointer to the second source matrix (9 floats).
void MulMat33(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL a,
  float const* MBASE_NOT_NULL b
);

/// ## Parameters
/// - `dst`: Pointer to the destination matrix (16 floats).
/// - `a`: Pointer to the first source matrix (16 floats).
/// - `b`: Pointer to the second source matrix (16 floats).
void MulMat44(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL a,
  float const* MBASE_NOT_NULL b
);

bool InvertMat33(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL src
);

bool InvertMat44(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL src
);

} // namespace mmath
