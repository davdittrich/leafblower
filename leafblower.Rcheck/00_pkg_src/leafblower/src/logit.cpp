#include "logit.hpp"
// LinkFn methods are all inline in the header; this TU is a placeholder
// for future non-inline implementations and ensures the TU is compiled.
// A static_assert on sizeof(LinkFn) would always pass (any struct has size >= 1)
// and is therefore not used here. Compilation of this TU itself serves as the
// "complete type" check.
namespace lbw {
} // namespace lbw
