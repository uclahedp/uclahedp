# uclahedp
UCLA HEDP Code

Analysis code for the ULCA HEDP group


# Installation (development mode)

The best way to install the package is in Anaconda's development mode. Clone the
repository onto your computer (ideally using GitHub Desktop) and then follow the 
instructions below to link the package to Anaconda. In this development mode, any changes
made to the local files (local edits or updates by pulling from github) will be
automatically reflected in the Anaconda package.

(Using Anaconda)
1. Open an Anaconda prompt (Windows) or a regular terminal prompt (Mac, Linux)
2. Run the command 'conda develop _' where _ is the package directory path. Eg.
'conda develop /Users/peter/Documents/GitHub/uclahedp'
3. Open Anaconda (close and restart if necessary) and test that the package has been
added to the path by importing something, eg. "from uclahedp import hdftools"
4. If necessary, packages can be uninstalled with 'conda develop -u _'