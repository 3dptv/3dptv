Prerequisities
---------------

1. Download and install ActiveTCL 8.4 ( not 8.5 or 8.6) e.g. 
http://downloads.activestate.com/ActiveTcl/releases/8.4.19.6/ActiveTcl8.4.19.6.295590-win32-ix86-threaded.exe

	- Install it to ``C:\Tcl``  (the path is important)



Download the pre-compiled binary for Windows 32- or 64-bit platform
---------------------


1. If you use 32bit  or 64 bit Windows platform you can use the binary distributions: 

	1. download the executable file from: https://github.com/3dptv/3dptv/downloads (Tcl/Tk 8.4 has to be installed), save it as ``C:\PTV\3dptv.exe``
	2. download the test folder https://github.com/downloads/3dptv/3dptv/test.zip, extract it to ``C:\PTV\test`` and test the software by running the ``C:\PTV\test\start.bat``

2. If you use another platform or you want to compile from the source, download the latest snapshot:
https://github.com/3dptv/3dptv/zipball/master


Compilation instructions from source on Windows platform
------------

1. If you want to compile from the source code on Windows, you recommend downloading the Visual C++  2010 Express http://www.microsoft.com/visualstudio/en-us/products/2010-editions/visual-cpp-express

	Important note for the first-time installation - RESTART YOUR WINDOWS after installing the ActiveTcl and Visual C++ 

2. Double-click the ``3dptv_vc2010.sln`` to open the 3DPTV package in Visual C++ 2010 Express. 

3. If all the paths are the same as ours, Build the solution and you should see the 3DPTV software up and running

If you use different paths it is important to:
	a. add ``\Tcl\Include`` to the list of included directories
	b. add ``\Tcl\tcl84.lib`` and ``\Tcl\tk84.lib`` to the additional libraries
	c. on Windows XP it is important to ignore ``LIBC.LIB`` and also to use ``MFC as Static Library``
	
	
	
Compilation from source on _nix platforms
-----------------------------------------

1. Read https://github.com/3dptv/3dptv/blob/master/install_nix.md


Multi-platform build using CMake
-----------------------------------------
1. Download and install CMake (version 2.8 or higher) 

	http://www.cmake.org/cmake/resources/software.html

2. Download and install ActiveTCL 8.4 for the desired platform
	
	http://www.activestate.com/activetcl/downloads

3. (linux systems) Download and install gcc and libtiff 

 	-- for Ubuntu run the following in a bash shell
	sudo apt-get install build-essential
	sudo apt-get install libtiff-dev 

4. create a build directory inside the main 3dptv directory (3dptv/build)


5. Run CMake and generate the project or make files 

	a. start the cmake-gui from the desktop shortcut or program list
	   - alternatively run the following in a bash shell
           cd

	b. browse for the source directory: 3dptv/

	c. browse for the build directory just created: 3dptv/build/

	d. click configure
	
	Note: if an error occurs saying that libtiff or ActiveTCL cannot be found
              click on the entry in the "value" column and browse for the libraries' "include" and "lib" directories
              Click configure again to locate the libraries. If any rows appear red, click configure once more.
	
	e. click generate to create the project files (Visual studio) or make files (linux)

6. Run Visual Studio and open the solution file in 3dptv/build/ or run make from within the build directory to compile the code.  
	--The executable will be generated and placed inside 3dptv/build/bin/