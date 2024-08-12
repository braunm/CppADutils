
# Notes for CppADutils setup and installation

- If you do not have Homebrew installed on your Mac, go to https://brew.sh to do that.

   Homebrew should install the Xcode Command Line tools for you.

- Install Eigen:

	    brew install eigen

(may not need to do that since we use RcppEigen anyway)

- Install CppAD:

	    brew install cppad

- In your .zprofile or .bash_profile (or wherever you set environment variables), add the following variable:

--  If you have an Intel Mac:

	export CPPAD_DIR=/usr/local/opt/cppad

-- If you have an Apple Silicon Mac:

	export CPPAD_DIR=/opt/homebrew/opt/cppad

Doing that will keep CppAD up-to-date in a fixed location, which is a lot easier than going through the standard CppAD installation instructions.

- Make sure R is at least version 4.1.  Update all of your R packages using update.packages()

- Download the lastest CppADutils:

        git clone https://github.com/braunm/CppADutils.git

make sure you are on the master branch

Launch R from the CppADutils package directory. and load the devtools package.  You may need to install the `pracma` package, along with some others.

	library(devtools)
	load_all(recompile=TRUE)
	test()

If all the tests pass, then install:

	install()

You will need to install so adshdlm can find CppADutils.

# Setting up adshdlm



To install adshdlm, update the latest from github:

	git clone https://github.com/braunm/adshdlm.git

There are a lot of branches, and I'm not sure which is the latest.  The most important update is that the src/Makevars file should look like:

    ifdef CPPAD_DIR
    PKG_CPPFLAGS= -I$(CPPAD_DIR)/include  -I../inst/include/ -DNDEBUG -DCPPAD_USE_CPLUSPLUS_2011 -Wno-unused-variable
    else
    PKG_CPPFLAGS=-I/opt/homebrew/opt/cppad/include  -I../inst/include/ -DNDEBUG -Wno-unused-variable
    endif
    PKG_CFLAGS=-O3 -Wno-unknown-pragmas
    CXX_STD = CXX11

If you are in a branch where src/Makevars is different, you can change it manually.  Or, go to the new branch in a command line and type

    git checkout --patch devel src/Makevars

That will copy the version of src/Makevars I updated in the devel branch and copy it to whichever branch you are in.

Then, do

	library(devtools)
	load_all(recompile=TRUE)
	test()

Hopefully, all the tests will pass.

Whenever you work, you can just do `load_all()` (don't always need to force recompiling).
