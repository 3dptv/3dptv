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




