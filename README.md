This document describes how to make existing ELiPS library working in Linux environment.
This is expected that it will work any 64bit Unix distribution Ububtu.
Autotools installation may vary for Linux.
Please keep in mind that it is still in developing phase. 
If found any bug related to installation, please infrom in `yuto_takahashi@s.okayama-u.ac.jp`


1. Follow the instructions to install `GMP` library. Latest vesion is ok.
2. Check if `autoconf` is installed in your environment. `autoconf --version`. You migh see someting like this . If it is not installed then follow point 3.
    autoconf (GNU Autoconf) 2.69
    Copyright (C) 2012 Free Software Foundation, Inc.
    License GPLv3+/Autoconf: GNU GPL version 3 or later
    <http://gnu.org/licenses/gpl.html>, <http://gnu.org/licenses/exceptions.html>
    This is free software: you are free to change and redistribute it.
    There is NO WARRANTY, to the extent permitted by law.
    Written by David J. MacKenzie and Akim Demaille.
    
3. Install `autoconf`  as follows 
    `sudo apt-get update`
    `sudo apt-get install autoconf`
4. Install `libtool` as follows 
    `sudo apt-get install libtool-bin`
5. git clone `https://github.com/ISecOkayamaUniv/elips_bls12_x64.git`
6. From terminal enter to `<elips_bls12_x64>` directory.
7. Run the following commands 
    `autoreconf -i`

The output will be almost as follows

    libtoolize: Consider adding 'AC_CONFIG_MACRO_DIRS([m4])' to configure.ac,
    libtoolize: and rerunning libtoolize and aclocal.
    libtoolize: Consider adding '-I m4' to ACLOCAL_AMFLAGS in Makefile.am.
    libtoolize: 'AC_PROG_RANLIB' is rendered obsolete by 'LT_INIT'
    configure.ac:7: warning: AM_INIT_AUTOMAKE: two- and three-arguments forms are deprecated.  For more info, see:
    configure.ac:7: http://www.gnu.org/software/automake/manual/automake.html#Modernize-AM_005fINIT_005fAUTOMAKE-invocation 
8. Next run `./configure` 
9. Then `make`
10. Finally `sudo make install`
10. To uninstall `sudo make uninstall` from the directory

Still there is no single header. Therefore please `/usr/local/include/ELiPS` directory to get the header declaration. 

If you face `cannot open shared object file: No such file or directory` while running then follow this steps:

1. Run from terminal 
   ` sudo ldconfig`
   
2. 
    `echo $LD_LIBRARY_PATH`

3. If the command of point 2 gives blank result then
    `LD_LIBRARY_PATH=/usr/local/lib`
4. Check if again of `echo $LD_LIBRARY_PATH`. If path is set then run again.

Retun to github account https://github.com/ISecOkayamaUniv 
