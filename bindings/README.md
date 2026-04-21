# Bindings

This directory is reserved for language packages that wrap the C++ core.

Planned subdirectories:
- `bindings/python/`
- `bindings/r/`

Recommended implementation approach:
- use `include/robustcause/capi.h` when you want the most stable FFI boundary
- use `include/robustcause/robustcause.hpp` when a richer C++ binding layer is worth the tighter coupling

Current work:
- `bindings/r/robustcause/`: self-contained R package wrapper around the C++ core
