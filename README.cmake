The CMake configuration for codec2 should be considered experimental at
this time but has been thouroughly tested on Fedora Linux and cross-compiling 
from linux to windows with mingw and has many advanages over the autotools
config. 

- Builds against system libraries (default).
- Has experimental NSIS packaing support for Windows (WIN32) targets. *nix
  systems should rely on 'make install' as the packages (RPM & DEB) created by
  CPack are questionable.

==========================
 Building and installing
==========================

    To test the cmake build make a directory anywhere underneath (or outside of)
the source directory.

Linux command line example:

$ cd /path/to/codec2
$ mkdir build_linux
$ cd build_linux
$ cmake ../ (defaults to /usr/local, use CMAKE_INSTALL_PREFIX to override)
(if no errors)
$ make
(as root)
$ make install
