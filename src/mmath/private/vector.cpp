// TU header --------------------------------------------
#include "mmath/public/vector.h"

// c++ headers ------------------------------------------
#include <cmath>

namespace mmath {

float Vector2::Length() const {
  return std::sqrt(this->Length2());
}

float Vector3::Length() const {
  return std::sqrt(this->Length2());
}

} // wentos::math
