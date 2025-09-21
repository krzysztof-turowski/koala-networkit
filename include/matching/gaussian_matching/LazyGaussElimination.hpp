#pragma once

#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>

#include <matching/gaussian_matching/utils.hpp>

namespace Koala {
class LazyGaussElimination {
 public:
  static std::vector<int>
  pivotElimination(MatZp &A, std::function<bool(int, int)> isCellAllowed);
  static std::vector<int> simpleElimination(MatZp &A, int k);
};
}  // namespace Koala
