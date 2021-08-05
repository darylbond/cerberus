#ifndef MFP_DEFS_H
#define MFP_DEFS_H

#ifdef AMREX_USE_EB
#define EB_OPTIONAL(...) __VA_ARGS__
#else
#define EB_OPTIONAL(...)
#endif

#endif // MFP_DEFS_H
