# ppm2_server_code

Fork of https://cggit.cc.lehigh.edu/biomembhub/ppm2_server_code with code modernization fixes.

It is the code behind the [Positioning of proteins in membranes (PPM) 
Web Server](https://opm.phar.umich.edu/ppm_server).

I am keeping the changes minimal, with a wrapper to workaround the special $pwd requirements. 

At the time of forking, the license was GPL v3.


## Testing

This should work

    make
    ./immers <test.inp >test.out

## References

* Lomize M.A., Pogozheva I,D, Joo H., Mosberg H.I., Lomize A.L. OPM database and PPM web server: resources for positioning of proteins in membranes. Nucleic Acids Res., 2012, 40 (Database issue):D370-376 
* Lomize AL, Todd SC, Pogozheva ID. (2022) Spatial arrangement of proteins in planar and curved membranes by PPM 3.0. Protein Sci. 31:209-220. 
