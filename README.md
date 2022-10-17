yzLibrary 
===================
A collection of useful structures and functions for easy development of demo projects. 

Author: Yizhong Zhang

<br/>

I originally write yzLibrary to ease the labor of developing small demos for computer graphics. I found it really disgusting to set dependencies every time I start a new project. It is even more disgusting to compile other people's code with a lot of third party libraries, even if these libraries are not necessary. So I write my own header only library. As more and more code added into yzLibrary, I reorganized the code structure and used [cmake](https://cmake.org/) to manage the code and its submodules. Now it is much easier to integrate into other projects. 

Currently, yzLibrary has 2 submodules, ./deps and 3rdParty/eigen. ./deps contains runtime dll for OpenGL and data necessary for example projects. 3rdParty/eigen is used in some functions. yzLibrary also works without these submodules, as you can use it just like a template only header library. 


## Usage

### Setup by cmake

It is recommended to use yzLibrary via cmake as submodule. In this way, all configs can be set automatically, and you can try the examples provided by yzLibrary very easily. 

Assume you have cloned yzLibrary as submodule in PRJ_ROOT/3rdParty/yzLibrary, add the following code into the project CMakeLists.txt file. By default, we set GLUT and GLEW enable in this code. You can manually change its value in cmake, or set other variable values accordingly.

<pre><code>
if(EXISTS ${CMAKE_SOURCE_DIR}/3rdParty/yzLibrary/yzLib)
    include(${CMAKE_SOURCE_DIR}/3rdParty/yzLibrary/yzLib/yzLib.cmake)
    yzLib_set_ENABLE(GLUT ON)
    yzLib_set_ENABLE(GLEW ON)
    
    include_directories(${yzLib_INCLUDE_DIRS})
    message(STATUS "yzLib included")
    yzLib_Print_Info()

    # ==================================================
    # the following 3 lines set dependencies for yzLibrary for runtime api, can be removed if you do not need them 
    include_directories(${CMAKE_SOURCE_DIR}/3rdParty/yzLibrary/deps/inc/)
    link_directories(${CMAKE_SOURCE_DIR}/3rdParty/yzLibrary/deps/lib/win64)
    set(BIN_PATH "${BIN_PATH};${CMAKE_SOURCE_DIR}/3rdParty/yzLibrary/deps/bin/win64")
    # ==================================================
else()
    message(FATAL_ERROR 
        "yzLib not exist, update git submodule to download. "
        "This can be achieved using the following command:\n "
        "git submodule update --recursive --init -- \"3rdParty/yzLibrary\"")
endif()
</code></pre>

After running this code, directory of yzLib header files will be included. The 3 lines between === will add dependencies of yzLibrary (OpenGL, etc.).

For more usage details, please refer the CMakeLists.txt in yzLibrary. This file will generate example projects for yzLibrary.


### Setup by copy code

If you prefer to use yzLibrary without cmake, you need to follow the following steps.

1. Copy the whole directory of yzLib to your include path.

2. Copy and rename yzLib_config.h.in as yzLib_config.h into your include path, then edit the file. For variables you need to set, change #cmakedefine into #define, otherwise remove the line.
