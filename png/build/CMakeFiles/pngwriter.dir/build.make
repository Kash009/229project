# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.0.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.0.0/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build

# Include any dependencies generated for this target.
include CMakeFiles/pngwriter.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/pngwriter.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/pngwriter.dir/flags.make

CMakeFiles/pngwriter.dir/src/pngwriter.cc.o: CMakeFiles/pngwriter.dir/flags.make
CMakeFiles/pngwriter.dir/src/pngwriter.cc.o: /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter/src/pngwriter.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/pngwriter.dir/src/pngwriter.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pngwriter.dir/src/pngwriter.cc.o -c /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter/src/pngwriter.cc

CMakeFiles/pngwriter.dir/src/pngwriter.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pngwriter.dir/src/pngwriter.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter/src/pngwriter.cc > CMakeFiles/pngwriter.dir/src/pngwriter.cc.i

CMakeFiles/pngwriter.dir/src/pngwriter.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pngwriter.dir/src/pngwriter.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter/src/pngwriter.cc -o CMakeFiles/pngwriter.dir/src/pngwriter.cc.s

CMakeFiles/pngwriter.dir/src/pngwriter.cc.o.requires:
.PHONY : CMakeFiles/pngwriter.dir/src/pngwriter.cc.o.requires

CMakeFiles/pngwriter.dir/src/pngwriter.cc.o.provides: CMakeFiles/pngwriter.dir/src/pngwriter.cc.o.requires
	$(MAKE) -f CMakeFiles/pngwriter.dir/build.make CMakeFiles/pngwriter.dir/src/pngwriter.cc.o.provides.build
.PHONY : CMakeFiles/pngwriter.dir/src/pngwriter.cc.o.provides

CMakeFiles/pngwriter.dir/src/pngwriter.cc.o.provides.build: CMakeFiles/pngwriter.dir/src/pngwriter.cc.o

# Object files for target pngwriter
pngwriter_OBJECTS = \
"CMakeFiles/pngwriter.dir/src/pngwriter.cc.o"

# External object files for target pngwriter
pngwriter_EXTERNAL_OBJECTS =

libpngwriter.dylib: CMakeFiles/pngwriter.dir/src/pngwriter.cc.o
libpngwriter.dylib: CMakeFiles/pngwriter.dir/build.make
libpngwriter.dylib: /usr/lib/libz.dylib
libpngwriter.dylib: /usr/local/lib/libpng.dylib
libpngwriter.dylib: /usr/lib/libz.dylib
libpngwriter.dylib: /usr/local/lib/libfreetype.dylib
libpngwriter.dylib: /usr/local/lib/libpng.dylib
libpngwriter.dylib: /usr/local/lib/libfreetype.dylib
libpngwriter.dylib: CMakeFiles/pngwriter.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library libpngwriter.dylib"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pngwriter.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/pngwriter.dir/build: libpngwriter.dylib
.PHONY : CMakeFiles/pngwriter.dir/build

CMakeFiles/pngwriter.dir/requires: CMakeFiles/pngwriter.dir/src/pngwriter.cc.o.requires
.PHONY : CMakeFiles/pngwriter.dir/requires

CMakeFiles/pngwriter.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/pngwriter.dir/cmake_clean.cmake
.PHONY : CMakeFiles/pngwriter.dir/clean

CMakeFiles/pngwriter.dir/depend:
	cd /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build/CMakeFiles/pngwriter.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/pngwriter.dir/depend

