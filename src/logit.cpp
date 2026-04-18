#include "logit.hpp"
// LinkFn methods are all inline in the header; this TU is a placeholder
// for future non-inline implementations and ensures the TU is compiled.
namespace lbw {
// Force instantiation to catch any template/inline errors at compile time
static_assert(sizeof(LinkFn) > 0, "LinkFn must be a complete type");
} // namespace lbw
