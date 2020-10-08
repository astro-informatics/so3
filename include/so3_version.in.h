#ifndef SO3_VERSION_H
#define SO3_VERSION_H
inline const char *so3_version_string() { return "@PROJECT_VERSION@"; }
inline const char *so3_info() {
  return "package:\n"
         "  name: SO3\n"
         "  description: Fast and accurate Wigner transforms\n"
         "  authors:\n"
         "      - Jason McEwen\n"
         "      - Martin Buttner\n"
         "      - Boris Leistedt\n"
         "  license: GPL-3\n"
         "  url: https://astro-informatics.github.io/so3\n"
         "  version: @PROJECT_VERSION@\n";
};
// clang-format off
inline int so3_version_major() { return @PROJECT_VERSION_MAJOR@; }
inline int so3_version_minor() { return @PROJECT_VERSION_MINOR@; }
inline int so3_version_patch() { return @PROJECT_VERSION_PATCH@; }
#define SO3_VERSION_MAJOR  @PROJECT_VERSION_MAJOR@
#define SO3_VERSION_MINOR  @PROJECT_VERSION_MINOR@
#define SO3_VERSION_PATCH  @PROJECT_VERSION_PATCH@
// clang-format on
#endif
