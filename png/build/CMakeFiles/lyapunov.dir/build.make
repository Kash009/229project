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
include CMakeFiles/lyapunov.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/lyapunov.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/lyapunov.dir/flags.make

CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o: CMakeFiles/lyapunov.dir/flags.make
CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o: /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter/examples/lyapunov.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o -c /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter/examples/lyapunov.cc

CMakeFiles/lyapunov.dir/examples/lyapunov.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lyapunov.dir/examples/lyapunov.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter/examples/lyapunov.cc > CMakeFiles/lyapunov.dir/examples/lyapunov.cc.i

CMakeFiles/lyapunov.dir/examples/lyapunov.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lyapunov.dir/examples/lyapunov.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter/examples/lyapunov.cc -o CMakeFiles/lyapunov.dir/examples/lyapunov.cc.s

CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o.requires:
.PHONY : CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o.requires

CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o.provides: CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o.requires
	$(MAKE) -f CMakeFiles/lyapunov.dir/build.make CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o.provides.build
.PHONY : CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o.provides

CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o.provides.build: CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o

# Object files for target lyapunov
lyapunov_OBJECTS = \
"CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o"

# External object files for target lyapunov
lyapunov_EXTERNAL_OBJECTS =

lyapunov: CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o
lyapunov: CMakeFiles/lyapunov.dir/build.make
lyapunov: libpngwriter.a
lyapunov: /usr/local/lib/libfreetype.dylib
lyapunov: /usr/lib/libz.dylib
lyapunov: /usr/local/lib/libpng.dylib
lyapunov: /usr/lib/libz.dylib
lyapunov: /usr/local/lib/libpng.dylib
lyapunov: /usr/local/lib/libfreetype.dylib
lyapunov: CMakeFiles/lyapunov.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable lyapunov"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lyapunov.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/lyapunov.dir/build: lyapunov
.PHONY : CMakeFiles/lyapunov.dir/build

CMakeFiles/lyapunov.dir/requires: CMakeFiles/lyapunov.dir/examples/lyapunov.cc.o.requires
.PHONY : CMakeFiles/lyapunov.dir/requires

CMakeFiles/lyapunov.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/lyapunov.dir/cmake_clean.cmake
.PHONY : CMakeFiles/lyapunov.dir/clean

CMakeFiles/lyapunov.dir/depend:
	cd /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build/CMakeFiles/lyapunov.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/lyapunov.dir/depend

