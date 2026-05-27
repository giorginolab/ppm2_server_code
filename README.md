# ppm2_server_code

Fork of https://cggit.cc.lehigh.edu/biomembhub/ppm2_server_code with code modernization fixes.

It is the code behind the [Positioning of proteins in membranes (PPM) 
Web Server](https://opm.phar.umich.edu/ppm_server).

I am keeping the changes minimal, with a wrapper to workaround the special $pwd requirements. 

At the time of forking, the license was GPL v3.


## Test

This should work on most platforms:

    make
    ./immers <test.inp >test.out

and produce the files `1gzmout.pdb` and `1rsyout.pdb`, plus the (hopefully expected) 

    warning Note: The following floating-point exceptions are signalling: IEEE_INVALID_FLAG IEEE_UNDERFLOW_FLAG

The outputs for an M1 MacOS with gfortran 15.2.0 are in `ref`.
They differ from the original reference ones (now in `ref_orig`) by fractions of an Angstrom.

## References

* Lomize M.A., Pogozheva I,D, Joo H., Mosberg H.I., Lomize A.L. OPM database and PPM web server: resources for positioning of proteins in membranes. Nucleic Acids Res., 2012, 40 (Database issue):D370-376 
* Lomize AL, Todd SC, Pogozheva ID. (2022) Spatial arrangement of proteins in planar and curved membranes by PPM 3.0. Protein Sci. 31:209-220. 
