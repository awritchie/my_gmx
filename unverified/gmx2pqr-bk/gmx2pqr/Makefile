# Generated automatically from Makefile.in by configure.
#
# This is a Gromacs 3.0 gmx2pqr makefile for your own utility programs.
#
# Copy this file to whatever directory you are using for your own
# software and add more targets like the gmx2pqr one below.
#
# If you are using gmake it is relatively straightforward to add
# an include based on environment variables (like previous Gromacs versions)
# to select compiler flags and stuff automatically, but below it is static:
#

# Variables set by the configuration script:
LIBS         = -lmd -lgmx -lxml2 -lfftw3f -lm   
LDFLAGS      = -L/usr/local/gromacs/lib -I/usr/local/lib  -framework Accelerate -std=c99
CFLAGS	     = -O3 -fomit-frame-pointer -finline-functions -Wall -Wno-unused -m64 -I/usr/include/libxml2 -I/usr/local/gromacs/include -I/usr/local/gromacs/include/gromacs -I.
CC           = gcc
LD           = $(CC) 

# The real make targets - note that most make programs support
# the shortcut $^ instead of listing all object files a second
# time, but we cannot count on it...

gmx2pqr:	gmx2pqr.c
		$(LD) $(LDFLAGS) -o $@ gmx2pqr.c $(LIBS) $(CFLAGS) 

clean :
	rm gmx2pqr
