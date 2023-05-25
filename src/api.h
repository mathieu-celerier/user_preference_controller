#pragma once

#if defined _WIN32 || defined __CYGWIN__
#  define UserPreferenceController_DLLIMPORT __declspec(dllimport)
#  define UserPreferenceController_DLLEXPORT __declspec(dllexport)
#  define UserPreferenceController_DLLLOCAL
#else
// On Linux, for GCC >= 4, tag symbols using GCC extension.
#  if __GNUC__ >= 4
#    define UserPreferenceController_DLLIMPORT __attribute__((visibility("default")))
#    define UserPreferenceController_DLLEXPORT __attribute__((visibility("default")))
#    define UserPreferenceController_DLLLOCAL __attribute__((visibility("hidden")))
#  else
// Otherwise (GCC < 4 or another compiler is used), export everything.
#    define UserPreferenceController_DLLIMPORT
#    define UserPreferenceController_DLLEXPORT
#    define UserPreferenceController_DLLLOCAL
#  endif // __GNUC__ >= 4
#endif // defined _WIN32 || defined __CYGWIN__

#ifdef UserPreferenceController_STATIC
// If one is using the library statically, get rid of
// extra information.
#  define UserPreferenceController_DLLAPI
#  define UserPreferenceController_LOCAL
#else
// Depending on whether one is building or using the
// library define DLLAPI to import or export.
#  ifdef UserPreferenceController_EXPORTS
#    define UserPreferenceController_DLLAPI UserPreferenceController_DLLEXPORT
#  else
#    define UserPreferenceController_DLLAPI UserPreferenceController_DLLIMPORT
#  endif // UserPreferenceController_EXPORTS
#  define UserPreferenceController_LOCAL UserPreferenceController_DLLLOCAL
#endif // UserPreferenceController_STATIC