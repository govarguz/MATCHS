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

# Utility rule file for ug.

# Include the progress variables for this target.
include CMakeFiles/ug.dir/progress.make

CMakeFiles/ug: ../doc/ug/conf.py
CMakeFiles/ug: _espressopp.so
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/u/gvargas/code/e++SCv1/buildV1/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "HTML ug documentation available at /u/gvargas/code/e++SCv1/doc/ug/_build/html"
	cd /u/gvargas/code/e++SCv1/doc/ug && /cobra/u/system/soft/SLE_12_SP3/packages/x86_64/cmake/3.10.2/bin/cmake -E env PYTHONPATH=/u/gvargas/code/e++SCv1/buildV1:/u/gvargas/code/e++SCv1/buildV1/contrib: /mpcdf/soft/SLE_12_SP3/packages/x86_64/anaconda/2.5.1.0/bin/sphinx-build -b html . _build/html

ug: CMakeFiles/ug
ug: CMakeFiles/ug.dir/build.make

.PHONY : ug

# Rule to build all files generated by this target.
CMakeFiles/ug.dir/build: ug

.PHONY : CMakeFiles/ug.dir/build

CMakeFiles/ug.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ug.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ug.dir/clean

CMakeFiles/ug.dir/depend:
	cd /u/gvargas/code/e++SCv1/buildV1 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /u/gvargas/code/e++SCv1 /u/gvargas/code/e++SCv1 /u/gvargas/code/e++SCv1/buildV1 /u/gvargas/code/e++SCv1/buildV1 /u/gvargas/code/e++SCv1/buildV1/CMakeFiles/ug.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ug.dir/depend

