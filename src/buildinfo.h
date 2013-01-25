//  GBM by Greg Ridgeway  Copyright (C) 2003
//  License:    GNU GPL (version 2 or later)

#ifndef BUILDINFO_H
#define BUILDINFO_H

    #undef ERROR
    #include <R.h>

    #define GBM_FAILED(hr) ((unsigned long)hr != 0)
    typedef unsigned long GBMRESULT;
    #define GBM_OK 0
    #define GBM_FAIL 1
    #define GBM_INVALIDARG 2
    #define GBM_OUTOFMEMORY 3
    #define GBM_INVALID_DATA 4
    #define GBM_NOTIMPL 5

    #define LEVELS_PER_CHUNK ((unsigned long) 1)

    typedef unsigned long ULONG;
    typedef char *PCHAR;

    // #define NOISY_DEBUG

#endif // BUILDINFO_H
