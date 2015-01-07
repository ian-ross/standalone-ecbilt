standalone-ecbilt
=================

Standalone version of ECBILT atmosphere model for testing

To start with, compile `inputfiles.f` and run the resulting
executable.  This program reads ASCII files and converts these to
binary input datafiles.

Next run make to make the SST forced atmosphere only executable.

Possible problems: compliler options in the makefile; problems when
the program attempts to write to a direct access file without
specifying the record; namelist format.

Go to the `run` directory and type `../ecbilt' to start a 10 year
integration from an arbitrary initial state.  Go to
`run/outputdata/atmos` to postprocess the data from this run.

Compare the results with the files that were produced at KNMI.

In need of help? Mail: selten@knmi.nl
