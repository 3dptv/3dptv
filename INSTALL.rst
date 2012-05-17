Prerequisities
---------------

1. Download and install ActiveTCL 8.4 ( not 8.5 or 8.6) e.g. 
http://downloads.activestate.com/ActiveTcl/releases/8.4.19.6/ActiveTcl8.4.19.6.295590-win32-ix86-threaded.exe

- Install it to C:\\Tcl  (the path is important) ``fixed space font``

2. Download Visual C++  2010 Express http://www.microsoft.com/visualstudio/en-us/products/2010-editions/visual-cpp-express

- Install it anywhere, e.g. C:\\Program Files\\VC2010

3. Download and unpack the newest version of 3DPTV software from our tarball on Github:
https://github.com/3dptv/3dptv/zipball/master

- Install it to C:\\PTV (the path is important)

4. Double-click the C:\\PTV\\3dptv_vc2010\\3dptv_vc2010.sln to open the 3DPTV package in Visual C++ 2010 Express. 

5 If all the paths are the same as ours, Build the solution and you should see the 3DPTV software up and running

If you use different paths it is important to:
	a. add \\Tcl\\Include to the list of included directories
	b. add \\Tcl\\tcl84.lib and \Tcl\tk84.lib to the additional libraries
	c. on Windows XP it is important to ignore LIBC.LIB and also to use MFC as Static Library




