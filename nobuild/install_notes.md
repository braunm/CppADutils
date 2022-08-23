
# Notes for CppADutils setup and installation

- If you do not have Homebrew installed on your Mac, go to https://brew.sh to do that.

   Homebrew should install the Xcode Command Line tools for you.

- Install Eigen:  brew install eigen  (may not need to do that since we use RcppEigen anyway)

- Install CppAD:  brew install cppad

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

Launch R from the CppADutils package directory. and load the devtools package.  You may need to install the pracma package, along with some others.

library(devtools)
load_all(recompile=TRUE)
test()

If all the tests pass, then install:

install()


To install adshdlm, update the latest from github.
