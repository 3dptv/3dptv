## If you want to use the source code and contribute to the project:

1. Fork this repository to your account, using the Fork button
2. create branch or fix the master branch if you feel confident enough
3. send us the pull request when you feel like ready to submit your branch or your master into the community branch


 
## Multi-platform build using CMake

1. Download and install CMake (version 2.8 or higher) 

  	http://www.cmake.org/cmake/resources/software.html

2. Download and install ActiveTCL 8.6 for the desired platform
	
		http://www.activestate.com/activetcl/downloads

3. (Linux systems) Download and install gcc and libtiff 

 	-- for Ubuntu run the following in a bash shell
 	
		sudo apt-get install build-essential
		sudo apt-get install libtiff-dev
	
on Mac OS X is recommended to install libtiff using Homebrew:
		brew install libtiff

4. create a build directory inside the main 3dptv directory (```3dptv/build```)


5. Run CMake and generate the project or make files 

	a. start the ```cmake-gui``` from the desktop shortcut or program list
	   - alternatively run the following in a bash shell  
			```cd```
	b. browse for the source directory: ```3dptv/```

	c. browse for the build directory just created: ```3dptv/build/```

	d. click ```configure``` (on Mac OS X choose the following options: UNIX makefiles, native compilers)
	
	Note: if an error occurs saying that libtiff or ActiveTCL cannot be found
              click on the entry in the "value" column and browse for the libraries' "include" and "lib" directories
              Click ```configure``` again to locate the libraries. If any rows appear red, click configure once more. 
	
	e. click generate to create the project files (Visual studio) or make files (linux)

6. Run Visual Studio and open the solution file in 3dptv/build/ or run make from within the build directory to compile the code.  
	--The executable will be generated and placed inside 3dptv/build/bin/
	
7. Run the software from the ```test``` folder with a link to ```ptv.tcl``` as an input:

		cd 3dptv/test
		../build/bin/3dptv ../ptv.tcl
