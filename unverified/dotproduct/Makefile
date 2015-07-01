# Generated automatically from Makefile.in by configure.
#
# This is a Gromacs 3.0 dotproduct makefile for your own utility programs.
#
# Copy this file to whatever directory you are using for your own
# software and add more targets like the dotproduct one below.
#
# If you are using gmake it is relatively straightforward to add
# an include based on environment variables (like previous Gromacs versions)
# to select compiler flags and stuff automatically, but below it is static:
#

# Variables set by the configuration script:
LIBS         = -lmd -lgmx -lxml2  -lnsl -lfftw3f -lm   -lSM -lICE -lX11   
LDFLAGS      = -L/home/densign/programs/gmx404/lib   
CFLAGS	     = -O3 -fomit-frame-pointer -finline-functions -Wall -Wno-unused -funroll-all-loops   -I/usr/include/libxml2 -I/home/densign/programs/gmx404/include -I/home/densign/programs/gmx404/include/gromacs	
CC           = cc
LD           = $(CC)

# The real make targets - note that most make programs support
# the shortcut $^ instead of listing all object files a second
# time, but we cannot count on it...

dotproduct:	dotproduct.o
		$(LD) $(LDFLAGS) -o $@ dotproduct.o $(LIBS)
