yzLibrary 
===================
A collection of useful structures and functions for easy development of demo projects. 

Author: Yizhong Zhang

<br/>

I originally write yzLibrary to ease the labor of developing small demos for computer graphics. I found it really disgusting to set dependencies every time I start a new project. It is even more disgusting to compile other people's code with a lot of third party libraries, even if these libraries are not necessary. So I write my own header only library. 

As the more and more code added into yzLibrary, I reorganized the code structure and used cmake to manage the resources. Now it is much easier to use and   

###Usage

It is recommended to use yzLibrary as submodule

<pre><code>
if(EXISTS ${CMAKE_SOURCE_DIR}/3rdParty/yzLibrary/yzLib)
	include(${CMAKE_SOURCE_DIR}/3rdParty/yzLibrary/yzLib/yzLib.cmake)
	yzLib_set_ENABLE(GLUT ON)
	yzLib_set_ENABLE(GLEW ON)

	include_directories(${yzLib_INCLUDE_DIRS})
	message(STATUS "yzLib included")
	yzLib_Print_Info()

	# yzLibrary deps
	if(EXISTS ${CMAKE_SOURCE_DIR}/3rdParty/yzLibrary/deps/inc)
		include_directories(${CMAKE_SOURCE_DIR}/3rdParty/yzLibrary/deps/inc/)
		link_directories(${CMAKE_SOURCE_DIR}/3rdParty/yzLibrary/deps/lib/win64)
		set(BIN_PATH "${BIN_PATH};${CMAKE_SOURCE_DIR}/3rdParty/yzLibrary/deps/bin/win64")
	else()
		message(FATAL_ERROR 
			"yzLibrary/deps not exist, update git submodule recursively to download. ")
	endif()
else()
	message(FATAL_ERROR 
		"yzLib not exist, update git submodule to download. "
		"This can be achieved using the following command:\n "
		"git submodule update --recursive --init -- \"3rdParty/yzLibrary\"")
endif()
</code></pre>

