# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /u/system/soft/SLE_12_SP3/packages/x86_64/cmake/3.10.2/bin/cmake

# The command to remove a file.
RM = /u/system/soft/SLE_12_SP3/packages/x86_64/cmake/3.10.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /u/gvargas/code/e++SCv1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /u/gvargas/code/e++SCv1/buildV1

# Utility rule file for symlink.

# Include the progress variables for this target.
include CMakeFiles/symlink.dir/progress.make

symlink: CMakeFiles/symlink.dir/build.make
	/cobra/u/system/soft/SLE_12_SP3/packages/x86_64/cmake/3.10.2/bin/cmake -E create_symlink /u/gvargas/code/e++SCv1/src espressopp
.PHONY : symlink

# Rule to build all files generated by this target.
CMakeFiles/symlink.dir/build: symlink

.PHONY : CMakeFiles/symlink.dir/build

CMakeFiles/symlink.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/symlink.dir/cmake_clean.cmake
.PHONY : CMakeFiles/symlink.dir/clean

CMakeFiles/symlink.dir/depend:
	cd /u/gvargas/code/e++SCv1/buildV1 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /u/gvargas/code/e++SCv1 /u/gvargas/code/e++SCv1 /u/gvargas/code/e++SCv1/buildV1 /u/gvargas/code/e++SCv1/buildV1 /u/gvargas/code/e++SCv1/buildV1/CMakeFiles/symlink.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/symlink.dir/depend

