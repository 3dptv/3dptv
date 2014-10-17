## This part of readme file is for the scanning PTV only:

1. Credits:

    @article{  
    year={2005},  
    issn={0723-4864},  
	journal={Experiments in Fluids},  
	volume={39},  
	number={5},  
	doi={10.1007/s00348-005-0031-7},  
	title={3D scanning particle tracking velocimetry},  
	url={http://dx.doi.org/10.1007/s00348-005-0031-7},  
	publisher={Springer-Verlag},  
	keywords={Scanning; Tracking; Lagrangian},  
	author={Hoyer, Klaus and Holzner, Markus and Lüthi, Beat and Guala, Michele and Liberzon, Alexander and Kinzelbach, Wolfgang},  
	pages={923-934},  
	language={English}  
	}


	@article{0957-0233-25-6-065301,  
	  author={D Krug and M Holzner and B Lüthi and M Wolf and A Tsinober and W Kinzelbach},  
	  title={A combined scanning PTV/LIF technique to simultaneously measure the full velocity gradient tensor and the 3D density field},  
	  journal={Measurement Science and Technology},  
	  volume={25},  
	  number={6},  
	  pages={065301},  
	  url={http://stacks.iop.org/0957-0233/25/i=6/a=065301},  
	  year={2014}  
	}

 2. Additional explanations

 This version of Scanning PTV is branched from the (non-scanning) regular 3D-PTV software in 2005 or so and therefore includes some bugs that were fixed in the OpenPTV version. On the other hand, this software was developed futher by the 
 ETH Zurich group over the years, with contributions in the source code of the calibration, sequencing and tracking and post-processing. Please consult the document put forward by Dominik Krug in the ```/docs''' directory


## If you want to use the source code and contribute to the project:

1. Fork this repository to your account, using the Fork button
2. create branch or fix the master branch if you feel confident enough
3. send us the pull request when you feel like ready to submit your branch or your master into the community branch


 
## Multi-platform build using CMake

1. Download and install CMake (version 2.8 or higher) 

  	http://www.cmake.org/cmake/resources/software.html

2. Download and install ActiveTCL 8.6 for the desired platform
	
	http://www.activestate.com/activetcl/downloads

3. Download and install libtiff (*nix and OSX systems) 

 	-- on Ubuntu run the following in a terminal:
 	
		sudo apt-get install build-essential
		sudo apt-get install libtiff-dev
	
	-- on Mac OS X is recommended to install libtiff using Homebrew:

		brew install libtiff

4. Create a build directory inside the main 3dptv directory (```3dptv/build```)


5. Run CMake and generate the project or make files 

	a. start the ```cmake-gui``` from the desktop shortcut or program list
		
		-- alternatively run the following from command line  

		cd 3dptv/build
		cmake-gui ../			

	b. browse for the source directory: ```3dptv/```

	c. browse for the build directory just created: ```3dptv/build/```

	d. click ```configure``` and select the build type and compiler
 
		--on Mac OS X choose the following options: UNIX makefiles, native compilers

		--on Windows choose Visual Studio 2010 or 2008 (32-bit supported currently)

		--on *nix systems choose Unix make files or alternatively use Eclipse CTD and unix make files
	
	Note: if an error occurs saying that libtiff or ActiveTCL cannot be found
              click on the entry in the "value" column and browse for the libraries' "include" directories and ".lib" files
              click ```configure``` again to locate the libraries. If any rows appear red, click configure until no rows are highlighted.
	
	e. click generate to create the project or make files

6. Run Visual Studio and open the solution file in ```3dptv/build/``` or run make from within the build directory to compile the code.  

	--The executable will be generated and placed inside ```3dptv/build/bin/```
	
7. Run the software from the ```test``` folder with a link to ```ptv.tcl``` as an input:

		cd 3dptv/test
		../build/bin/3dptv ../ptv.tcl
