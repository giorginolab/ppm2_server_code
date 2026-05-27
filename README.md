# ppm2_server_code

Fork of https://cggit.cc.lehigh.edu/biomembhub/ppm2_server_code with code modernization fixes.

It is the code behind the [Positioning of proteins in membranes (PPM) 
Web Server](https://opm.phar.umich.edu/ppm_server).

I am keeping the changes minimal, with a wrapper to work around `immers`
expecting to run from a particular working directory layout.


## Test

This should work on most platforms:

    make
    ./immers <test.inp >test.out

and produce the files `1gzmout.pdb` and `1rsyout.pdb`, plus the (hopefully expected) warning 

    Note: The following floating-point exceptions are signalling: IEEE_INVALID_FLAG IEEE_UNDERFLOW_FLAG

The outputs for an M1 MacOS with gfortran 15.2.0 are in `ref`.
They differ from the original reference ones (now in `ref_orig`) by fractions of an Angstrom.

## Known numerical sensitivity

This legacy Fortran code is sensitive to compiler optimization. On Apple
Silicon with gfortran 15.2.0, comparing `-O` and `-O2` builds changes the
rounded `test.out` values for `1rsy.pdb` and changes several final `1gzm.pdb`
TM segment boundaries and energy/thickness summaries.

The origin of the `1gzm.pdb` difference is the profile table calculation in
`profile.f`, specifically `profile_point`. Mixed-optimization builds show that
compiling only the profile table generation at `-O` while leaving the rest of
the program at `-O2` restores the `-O` result. Doing the same for `opm.f`,
`locate.f`, `min.f`, `deftm.f`, `watface.f`, `tilting.f`, or only `ener_at`
does not.

Those small profile-table differences feed into the membrane-thickness scan in
`optim` in `opm.f`. Around the selected minimum the adjusted energy surface is
very flat. In the observed run:

    -O:  dmembr 31.0 -> -76.75543, 31.4 -> -76.74245
    -O2: dmembr 31.0 -> -76.73165, 31.4 -> -76.74550

So the selected scan point flips by one step. The final report prints
`thickn = dmin + 0.8`, which gives `31.8` for the `-O` result and `32.2` for
the `-O2` result.

There are broader numerical hazards in the codebase. Many loops use real-valued
Fortran `DO` variables, which modern gfortran reports as a deleted feature.
The normal test run also leaves `IEEE_INVALID_FLAG` set. With floating-point
invalid traps enabled, one trigger is `acos(co)` in `watface.f`, where roundoff
can place `co` slightly outside the valid `[-1, 1]` range. This is separate
from the `test.out` versus `testO2.out` origin above, but it explains part of
the floating-point warning already noted in the test output.

## Wrapper

Use `run_immers_ppm2.py` to run the bundled `immers` executable on a single
PDB file without copying the input into `ppm2_server_code/` by hand:

    python run_immers_ppm2.py path/to/protein.pdb out

The wrapper creates a temporary working directory for each run, copies the input
PDB there as the fixed short name `p001.pdb`, links or copies `immers` and
`res.lib`, writes the PPM2 input file expected by `immers`, and runs:

    ./immers < input.inp

After the run, it copies new output files back next to the original PDB. Outputs
whose names contain the temporary stem `p001` are renamed to use the original
PDB stem, so an input like `c0.pdb` receives outputs such as `c0out.pdb`.
Table outputs `datapar1` and `datasub1` are copied back as `<stem>.datapar1`
and `<stem>.datasub1`, and their internal references are rewritten from
`p001.pdb` back to the original PDB name.

The generated `<stem>.inp` and `<stem>.opm2.stdout` files are also kept next to
the input PDB for reproducibility and debugging. Pass `--debug` to keep the
temporary directory and print its path.

At the time of forking, the license was GPL v3.



## References

* Lomize M.A., Pogozheva I,D, Joo H., Mosberg H.I., Lomize A.L. OPM database and PPM web server: resources for positioning of proteins in membranes. Nucleic Acids Res., 2012, 40 (Database issue):D370-376 
* Lomize AL, Todd SC, Pogozheva ID. (2022) Spatial arrangement of proteins in planar and curved membranes by PPM 3.0. Protein Sci. 31:209-220. 
