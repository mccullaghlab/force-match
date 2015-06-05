# Force_Match
Force Matching Python script and config file.

To use this program, you should only edit "forcematchFiles.config". The Python script runs
based on how the .config file is filled out.

Place any .psf, .pdb, .dcd, and .force.dcd files of interest in the program's directory. In the .config file,
name the specific set of .psf, .dcd, and .force.dcd files which are to be analyzed the next time the program is run, and save the file.

For example:

    PSF FILE:
    nacl2_wwat.psf
    
    PDB FILE:
    NONE
    
    FORCE DCD FILES:
    nacl2_wwat.md.run2.out.force.dcd
    
    COORD DCD FILES:
    nacl2_wwat.md.run2.out.dcd
    
    DEBUG MODE: OFF
    
    END CONFIG

As it says in the config file, only edit the file names and config settings. Do not edit the heading
lines, such as "PSF FILE:".
