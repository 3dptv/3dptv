Prerequisities
---------------

1. Download and install ActiveTCL 8.4 ( not 8.5 or 8.6) e.g. 
http://downloads.activestate.com/ActiveTcl/releases/8.4.19.6/ActiveTcl8.4.19.6.295590-win32-ix86-threaded.exe

- Install it to ``C:\Tcl``  (the path is important)


2. If you use 32bit Windows platform you can try to use the binary distribution:

	1. download https://github.com/downloads/3dptv/3dptv/3dptv.exe (Tcl/Tk 8.4 has to be installed)
	2. download the software package https://github.com/3dptv/3dptv/zipball/master
	3. uncompress the software package e.g. in C:\PTV and run the ``C:\PTV\test\start.bat``




3. If you need to compile from source, download Visual C++  2010 Express http://www.microsoft.com/visualstudio/en-us/products/2010-editions/visual-cpp-express

- Install it anywhere, e.g. ``C:\Program Files\VC2010``


Important note for the first-time installation - RESTART YOUR WINDOWS after installing the ActiveTcl and Visual C++ 


Download the software
---------------------

1. If you wish to compile from source - download and unpack the newest version of 3DPTV software from our tarball on Github:
https://github.com/3dptv/3dptv/zipball/master

2. If you want to try the pre-compiled binary executables visit https://github.com/3dptv/3dptv/downloads for the 
Windows 32bit XP or Windows 7 64bit (with Tcl/Tk 32bit libraries) versions. 

- Expand the compressed file  into ``C:\PTV`` folder (the path is important it is prescribed by the Visual C++ project file)

Compilation from source on Windows platform
------------

1. Double-click the ``C:\PTV\3dptv_vc2010\3dptv_vc2010.sln`` to open the 3DPTV package in Visual C++ 2010 Express. 

2. If all the paths are the same as ours, Build the solution and you should see the 3DPTV software up and running

If you use different paths it is important to:
	a. add ``\Tcl\Include`` to the list of included directories
	b. add ``\Tcl\tcl84.lib`` and ``\Tcl\tk84.lib`` to the additional libraries
	c. on Windows XP it is important to ignore ``LIBC.LIB`` and also to use ``MFC as Static Library``
	
	
	
Compilation from source on _nix platforms
-----------------------------------------

1. Read https://github.com/3dptv/3dptv/blob/master/install_nix.rst




