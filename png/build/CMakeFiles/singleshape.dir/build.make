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
include CMakeFiles/singleshape.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/singleshape.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/singleshape.dir/flags.make

CMakeFiles/singleshape.dir/examples/singleshape.cc.o: CMakeFiles/singleshape.dir/flags.make
CMakeFiles/singleshape.dir/examples/singleshape.cc.o: /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter/examples/singleshape.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/singleshape.dir/examples/singleshape.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/singleshape.dir/examples/singleshape.cc.o -c /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter/examples/singleshape.cc

CMakeFiles/singleshape.dir/examples/singleshape.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/singleshape.dir/examples/singleshape.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter/examples/singleshape.cc > CMakeFiles/singleshape.dir/examples/singleshape.cc.i

CMakeFiles/singleshape.dir/examples/singleshape.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/singleshape.dir/examples/singleshape.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter/examples/singleshape.cc -o CMakeFiles/singleshape.dir/examples/singleshape.cc.s

CMakeFiles/singleshape.dir/examples/singleshape.cc.o.requires:
.PHONY : CMakeFiles/singleshape.dir/examples/singleshape.cc.o.requires

CMakeFiles/singleshape.dir/examples/singleshape.cc.o.provides: CMakeFiles/singleshape.dir/examples/singleshape.cc.o.requires
	$(MAKE) -f CMakeFiles/singleshape.dir/build.make CMakeFiles/singleshape.dir/examples/singleshape.cc.o.provides.build
.PHONY : CMakeFiles/singleshape.dir/examples/singleshape.cc.o.provides

CMakeFiles/singleshape.dir/examples/singleshape.cc.o.provides.build: CMakeFiles/singleshape.dir/examples/singleshape.cc.o

# Object files for target singleshape
singleshape_OBJECTS = \
"CMakeFiles/singleshape.dir/examples/singleshape.cc.o"

# External object files for target singleshape
singleshape_EXTERNAL_OBJECTS =

singleshape: CMakeFiles/singleshape.dir/examples/singleshape.cc.o
singleshape: CMakeFiles/singleshape.dir/build.make
singleshape: libpngwriter.a
singleshape: /usr/local/lib/libfreetype.dylib
singleshape: /usr/lib/libz.dylib
singleshape: /usr/local/lib/libpng.dylib
singleshape: /usr/lib/libz.dylib
singleshape: /usr/local/lib/libpng.dylib
singleshape: /usr/local/lib/libfreetype.dylib
singleshape: CMakeFiles/singleshape.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable singleshape"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/singleshape.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/singleshape.dir/build: singleshape
.PHONY : CMakeFiles/singleshape.dir/build

CMakeFiles/singleshape.dir/requires: CMakeFiles/singleshape.dir/examples/singleshape.cc.o.requires
.PHONY : CMakeFiles/singleshape.dir/requires

CMakeFiles/singleshape.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/singleshape.dir/cmake_clean.cmake
.PHONY : CMakeFiles/singleshape.dir/clean

CMakeFiles/singleshape.dir/depend:
	cd /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/pngwriter /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build /Users/winnielin/Schoolwork/2015/Fall/CS229/Project/png/build/CMakeFiles/singleshape.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/singleshape.dir/depend

