# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /opt/local/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/LNLB/mpi-projects/FD2D

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/LNLB/mpi-projects/FD2D

# Include any dependencies generated for this target.
include fdmodels/CMakeFiles/models.dir/depend.make

# Include the progress variables for this target.
include fdmodels/CMakeFiles/models.dir/progress.make

# Include the compile flags for this target's objects.
include fdmodels/CMakeFiles/models.dir/flags.make

# Object files for target models
models_OBJECTS =

# External object files for target models
models_EXTERNAL_OBJECTS =

fdmodels/libmodels.a: fdmodels/CMakeFiles/models.dir/build.make
fdmodels/libmodels.a: fdmodels/CMakeFiles/models.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libmodels.a"
	cd /Users/LNLB/mpi-projects/FD2D/fdmodels && $(CMAKE_COMMAND) -P CMakeFiles/models.dir/cmake_clean_target.cmake
	cd /Users/LNLB/mpi-projects/FD2D/fdmodels && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/models.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
fdmodels/CMakeFiles/models.dir/build: fdmodels/libmodels.a
.PHONY : fdmodels/CMakeFiles/models.dir/build

fdmodels/CMakeFiles/models.dir/requires:
.PHONY : fdmodels/CMakeFiles/models.dir/requires

fdmodels/CMakeFiles/models.dir/clean:
	cd /Users/LNLB/mpi-projects/FD2D/fdmodels && $(CMAKE_COMMAND) -P CMakeFiles/models.dir/cmake_clean.cmake
.PHONY : fdmodels/CMakeFiles/models.dir/clean

fdmodels/CMakeFiles/models.dir/depend:
	cd /Users/LNLB/mpi-projects/FD2D && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/LNLB/mpi-projects/FD2D /Users/LNLB/mpi-projects/FD2D/fdmodels /Users/LNLB/mpi-projects/FD2D /Users/LNLB/mpi-projects/FD2D/fdmodels /Users/LNLB/mpi-projects/FD2D/fdmodels/CMakeFiles/models.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : fdmodels/CMakeFiles/models.dir/depend

