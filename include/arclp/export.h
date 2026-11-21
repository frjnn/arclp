/* Symbol visibility of the public API.  Shared builds compile everything
 * else hidden; ARCLP_SHARED is set for consumers by the CMake target. */
#ifndef ARCLP_EXPORT_H
#define ARCLP_EXPORT_H

#if defined(_WIN32) && defined(ARCLP_SHARED)
#if defined(ARCLP_BUILDING)
#define ARCLP_API __declspec(dllexport)
#else
#define ARCLP_API __declspec(dllimport)
#endif
#elif defined(__GNUC__) || defined(__clang__)
#define ARCLP_API __attribute__((visibility("default")))
#else
#define ARCLP_API
#endif

#endif
