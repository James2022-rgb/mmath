// TU header --------------------------------------------
#include "mmath/public/matrix.h"

// c++ headers ------------------------------------------
#include <cstring>
#include <cmath>

#include <memory>
#include <utility>

namespace mmath {

//
// Matrix33
//

Matrix33::Matrix33(float const* MBASE_NOT_NULL v) {
  memcpy(data.data(), v, 9 * sizeof(float));
}

Matrix33::Matrix33(
  float const* MBASE_NOT_NULL col0,
  float const* MBASE_NOT_NULL col1,
  float const* MBASE_NOT_NULL col2
) {
  memcpy(data.data() + 0,  col0, 3 * sizeof(float)); // Column 0
  memcpy(data.data() + 3,  col1, 3 * sizeof(float)); // Column 1
  memcpy(data.data() + 6,  col2, 3 * sizeof(float)); // Column 2
}

mmath::Matrix33 Matrix33::Transpose() const {
  Matrix33 result = *this;
  TransposeMat33Inplace(result.data.data());
  return result;
}
bool Matrix33::Invert(Matrix33& out_result) const {
  return InvertMat33(out_result.data.data(), this->data.data());
}

Matrix33 Matrix33::operator*(float scalar) const {
  Matrix33 result;
  for (size_t i = 0; i < 9; ++i) {
    result.data[i] = this->data[i] * scalar;
  }
  return result;
}

Matrix33 Matrix33::operator*(Matrix33 const& other) const {
  Matrix33 result;
  MulMat33(result.data.data(), this->data.data(), other.data.data());
  return result;
}

Vector3 operator*(Matrix33 const& mat, Vector3 const& vec) {
  Vector3 result;
  MulMat33Vec3(result.Data(), mat.data.data(), vec.Data());
  return result;
}

//
// Matrix44
//

Matrix44::Matrix44(float const* MBASE_NOT_NULL v) {
  memcpy(data.data(), v, 16 * sizeof(float));
}

Matrix44::Matrix44(
  float const* MBASE_NOT_NULL col0,
  float const* MBASE_NOT_NULL col1,
  float const* MBASE_NOT_NULL col2,
  float const* MBASE_NOT_NULL col3
) {
  memcpy(data.data() + 0,  col0, 4 * sizeof(float)); // Column 0
  memcpy(data.data() + 4,  col1, 4 * sizeof(float)); // Column 1
  memcpy(data.data() + 8,  col2, 4 * sizeof(float)); // Column 2
  memcpy(data.data() + 12, col3, 4 * sizeof(float)); // Column 3
}

Matrix44 Matrix44::Transpose() const {
  Matrix44 result = *this;
  TransposeMat44Inplace(result.data.data());
  return result;
}

bool Matrix44::Invert(Matrix44& out_result) const {
  return InvertMat44(out_result.data.data(), this->data.data());
}

Matrix44 Matrix44::operator*(float scalar) const {
  Matrix44 result;
  for (size_t i = 0; i < 16; ++i) {
    result.data[i] = this->data[i] * scalar;
  }
  return result;
}

Matrix44 Matrix44::operator*(Matrix44 const& other) const {
  Matrix44 result;
  MulMat44(result.data.data(), this->data.data(), other.data.data());
  return result;
}

Vector4 operator*(Matrix44 const& mat, Vector4 const& vec) {
  Vector4 result;
  MulMat44Vec4(result.Data(), mat.data.data(), vec.Data());
  return result;
}

void IdentityMat33(float* MBASE_NOT_NULL dst) {
  dst[0] = 1.0f; dst[1] = 0.0f; dst[2] = 0.0f;
  dst[3] = 0.0f; dst[4] = 1.0f; dst[5] = 0.0f;
  dst[6] = 0.0f; dst[7] = 0.0f; dst[8] = 1.0f;
}

void IdentityMat44(float* MBASE_NOT_NULL dst) {
  dst[0] = 1.0f; dst[1] = 0.0f; dst[2] = 0.0f; dst[3] = 0.0f;
  dst[4] = 0.0f; dst[5] = 1.0f; dst[6] = 0.0f; dst[7] = 0.0f;
  dst[8] = 0.0f; dst[9] = 0.0f; dst[10] = 1.0f; dst[11] = 0.0f;
  dst[12] = 0.0f; dst[13] = 0.0f; dst[14] = 0.0f; dst[15] = 1.0f;
}

void TransposeMat44(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL src
) {
  dst[0] = src[0]; // Diagonal
  dst[1] = src[4];
  dst[2] = src[8];
  dst[3] = src[12];
  dst[4] = src[1];
  dst[5] = src[5]; // Diagonal
  dst[6] = src[9];
  dst[7] = src[13];
  dst[8] = src[2];
  dst[9] = src[6];
  dst[10] = src[10]; // Diagonal
  dst[11] = src[14];
  dst[12] = src[3];
  dst[13] = src[7];
  dst[14] = src[11];
  dst[15] = src[15]; // Diagonal
}

void TransposeMat33Inplace(float* MBASE_NOT_NULL dst) {
  std::swap(dst[1], dst[3]);
  std::swap(dst[2], dst[6]);
  std::swap(dst[5], dst[7]);
}


void TransposeMat44Inplace(float* MBASE_NOT_NULL dst) {
  std::swap(dst[1], dst[4]);
  std::swap(dst[2], dst[8]);
  std::swap(dst[3], dst[12]);
  std::swap(dst[6], dst[9]);
  std::swap(dst[7], dst[13]);
  std::swap(dst[11], dst[14]);
}

void MulMat33Vec3(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL mat,
  float const* MBASE_NOT_NULL vec
) {
  dst[0] = mat[0] * vec[0] + mat[3] * vec[1] + mat[6] * vec[2];
  dst[1] = mat[1] * vec[0] + mat[4] * vec[1] + mat[7] * vec[2];
  dst[2] = mat[2] * vec[0] + mat[5] * vec[1] + mat[8] * vec[2];
}

void MulMat44Vec4(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL mat,
  float const* MBASE_NOT_NULL vec
) {
  dst[0] = mat[0] * vec[0] + mat[4] * vec[1] + mat[8]  * vec[2] + mat[12] * vec[3];
  dst[1] = mat[1] * vec[0] + mat[5] * vec[1] + mat[9]  * vec[2] + mat[13] * vec[3];
  dst[2] = mat[2] * vec[0] + mat[6] * vec[1] + mat[10] * vec[2] + mat[14] * vec[3];
  dst[3] = mat[3] * vec[0] + mat[7] * vec[1] + mat[11] * vec[2] + mat[15] * vec[3];
}

void MulMat33(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL a,
  float const* MBASE_NOT_NULL b
) {
  // Column 0
  dst[0] = a[0] * b[0] + a[3] * b[1] + a[6] * b[2];
  dst[1] = a[1] * b[0] + a[4] * b[1] + a[7] * b[2];
  dst[2] = a[2] * b[0] + a[5] * b[1] + a[8] * b[2];

  // Column 1
  dst[3] = a[0] * b[3] + a[3] * b[4] + a[6] * b[5];
  dst[4] = a[1] * b[3] + a[4] * b[4] + a[7] * b[5];
  dst[5] = a[2] * b[3] + a[5] * b[4] + a[8] * b[5];

  // Column 2
  dst[6] = a[0] * b[6] + a[3] * b[7] + a[6] * b[8];
  dst[7] = a[1] * b[6] + a[4] * b[7] + a[7] * b[8];
  dst[8] = a[2] * b[6] + a[5] * b[7] + a[8] * b[8];
}

void MulMat44(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL a,
  float const* MBASE_NOT_NULL b
) {
  // Row 0
  dst[0]  = a[0] * b[0]  + a[4] * b[1]  + a[8] * b[2]  + a[12] * b[3];
  dst[4]  = a[0] * b[4]  + a[4] * b[5]  + a[8] * b[6]  + a[12] * b[7];
  dst[8]  = a[0] * b[8]  + a[4] * b[9]  + a[8] * b[10] + a[12] * b[11];
  dst[12] = a[0] * b[12] + a[4] * b[13] + a[8] * b[14] + a[12] * b[15];

  // Row 1
  dst[1]  = a[1] * b[0]  + a[5] * b[1]  + a[9] * b[2]  + a[13] * b[3];
  dst[5]  = a[1] * b[4]  + a[5] * b[5]  + a[9] * b[6]  + a[13] * b[7];
  dst[9]  = a[1] * b[8]  + a[5] * b[9]  + a[9] * b[10] + a[13] * b[11];
  dst[13] = a[1] * b[12] + a[5] * b[13] + a[9] * b[14] +a[13]  * b[15];

  // Row 2
  dst[2]  = a[2] * b[0]  + a[6] * b[1]  + a[10] * b[2]  + a[14] * b[3];
  dst[6]  = a[2] * b[4]  + a[6] * b[5]  + a[10] * b[6]  + a[14] * b[7];
  dst[10] = a[2] * b[8]  + a[6] * b[9]  + a[10] * b[10] + a[14] * b[11];
  dst[14] = a[2] * b[12] + a[6] * b[13] + a[10] * b[14] + a[14] * b[15];

  // Row 3
  dst[3]  = a[3] * b[0]  + a[7] * b[1]  + a[11] * b[2]  + a[15] * b[3];
  dst[7]  = a[3] * b[4]  + a[7] * b[5]  + a[11] * b[6]  + a[15] * b[7];
  dst[11] = a[3] * b[8]  + a[7] * b[9]  + a[11] * b[10] + a[15] * b[11];
  dst[15] = a[3] * b[12] + a[7] * b[13] + a[11] * b[14] + a[15] * b[15];
}

bool InvertMat33(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL src
) {
  // Prevent accidental in-place usage
  float tmp[9];
  if (dst == src) {
    std::memcpy(tmp, src, sizeof(tmp));
    src = tmp;
  }

  // Precompute minors (2x2 determinants)
  float m00 = src[4] * src[8] - src[5] * src[7];
  float m01 = src[3] * src[8] - src[5] * src[6];
  float m02 = src[3] * src[7] - src[4] * src[6];
  float m10 = src[1] * src[8] - src[2] * src[7];
  float m11 = src[0] * src[8] - src[2] * src[6];
  float m12 = src[0] * src[7] - src[1] * src[6];
  float m20 = src[1] * src[5] - src[2] * src[4];
  float m21 = src[0] * src[5] - src[2] * src[3];
  float m22 = src[0] * src[4] - src[1] * src[3];

  float const det =
    src[0] * m00 -
    src[3] * m10 +
    src[6] * m20;

  constexpr float eps = 1e-8f;
  if (std::fabs(det) < eps) {
    return false;
  }

  float inv_det = 1.0f / det;

  // Cofactor matrix transposed directly into column-major inverse
  dst[0] =  m00 * inv_det;
  dst[1] = -m10 * inv_det;
  dst[2] =  m20 * inv_det;

  dst[3] = -m01 * inv_det;
  dst[4] =  m11 * inv_det;
  dst[5] = -m21 * inv_det;

  dst[6] =  m02 * inv_det;
  dst[7] = -m12 * inv_det;
  dst[8] =  m22 * inv_det;

  return true;
}

bool InvertMat44(
  float* MBASE_NOT_NULL dst,
  float const* MBASE_NOT_NULL src
) {
  // Convert column-major src to row-major augmented [A | I]
  float m[4 * 8]; // 4 rows, 8 columns
  for (int r = 0; r < 4; ++r) {
    for (int c = 0; c < 4; ++c) {
      m[r * 8 + c] = src[c * 4 + r]; // src(col-major) -> row-major
    }
    for (int c = 0; c < 4; ++c) {
      m[r * 8 + 4 + c] = (r == c) ? 1.0f : 0.0f;
    }
  }

  // Gauss-Jordan with partial pivoting
  for (int col = 0; col < 4; ++col) {
    // Find pivot row
    int pivotRow = col;
    float maxAbs = std::fabs(m[pivotRow * 8 + col]);
    for (int r = col + 1; r < 4; ++r) {
      float v = std::fabs(m[r * 8 + col]);
      if (v > maxAbs) {
        maxAbs = v;
        pivotRow = r;
      }
    }
    if (maxAbs < 1e-8f) {
      return false; // Singular (or near-singular)
    }
    // Swap pivot row
    if (pivotRow != col) {
      for (int k = 0; k < 8; ++k) {
        std::swap(m[col * 8 + k], m[pivotRow * 8 + k]);
      }
    }
    // Normalize pivot row
    float pivot = m[col * 8 + col];
    float invPivot = 1.0f / pivot;
    for (int k = 0; k < 8; ++k) {
      m[col * 8 + k] *= invPivot;
    }
    // Eliminate other rows
    for (int r = 0; r < 4; ++r) {
      if (r == col) continue;
      float factor = m[r * 8 + col];
      if (factor != 0.0f) {
        for (int k = 0; k < 8; ++k) {
          m[r * 8 + k] -= factor * m[col * 8 + k];
        }
      }
    }
  }

  // Extract inverse: convert back to column-major
  for (int r = 0; r < 4; ++r) {
    for (int c = 0; c < 4; ++c) {
      dst[c * 4 + r] = m[r * 8 + 4 + c];
    }
  }
  return true;
}

} // namespace mmath
