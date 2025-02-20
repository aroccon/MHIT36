# MHIT36_cuDecomp

Tentative porting of MHIT36 to multi GPU using cuDecomp.

Log of changes/status of the porting
- 06/02/25: MPI init in MHIT36 
- 07/02/25: Build of cuDecomp ok, present in cuDecomp-main/Build
- 11/02/25: Up to cuFFT plan, adding comments step by step 
- 12/02/25: Compile and run on Leonardo + makefile for Leo
- 13/02/25: Adaption to work with only one grid descriptor and pencils along y and z only (0,1,1). Everyhing seems fine
- 14/02/25: Pencil update on y and z look also fine and working, all basic elements seems to be ready. cuDecomp updated to (#52), version 0.4.2, Tested on Leonardo up to 16 nodes (with the makefile with -cuda + -Wl, relax). The makefile issue is present only on Leonardo where cuda is not present in the same folder.
- 15/02/25: General structure of the code inializated and division in substeps, rhsp and p as input and output for Poisson.
- 17/02/25: Cleaning and use of module.f90 and readinput.f90, more variable included in the modules. Temporal loop created. Soon merging with main repository?
- 18/02/25: Fixed issue on makefile and forced linking (working also in Leonardo), flow field initialization looks of for different PR and PC; MPI output in parallel; tested different grid resolutions (up to 1024^3).
- 19/02/25: Phase-field initialization and MPI I/O of phase variables seem fine. Implementation of the projection step (convective, diffusive and forcing); halo updates of ustar, vstar and wstar also implemented. Issue with 1536^3 and 2048^3 on 8 nodes (32 GPUs).
poisson.f90 has the same issue (which is the one provided by Nvidia).
- 20/02/25: Problem on large grid (1536^3 and 2048^3) has been fixed (Thank you Josh), there was an integer overflow in the normalization; projection step implemented; first run on Local machine of the full NS solver. Still something off in the solver; blows up after a few iterations.
Strange behaviof of the convective terms, introduced internal update if pr=1 or pc=1. Now using dx=lx/(nx-1) the convective terms are fine.

# Multi-GPU version status

- Poisson solver (transposition + halo update) ✅
- Poisson solver validation (periodic solutions) ✅
- Read input files ✅
- Skeleton of the code  ✅
- Halo updates test with CUDA ✅
- Poisson solver scaling ✅
- Halo updates test with host_data use_device ✅
- Flow field initialization ✅
- Phase-field initialization ✅
- Projection step implemented ✅
- Validation of projection step ✅ (implemented, not validated)
- Correction step ✅ (implemented, not validated)
- Forcing ✅
- HIT validation ❌
- Drop oscillation validation ❌
- Full code scaling ❌
- MPI writing (no halo)  ✅
- Serial reading (to avoid issue with Leonardo) ❌
- Courant number check (MPI reduction) ❌


# Run the code

- Compile first the cuDecomp library using *_lib.sh, the resulting modules and library will be located in cuDecomp/build/lib and cuDecomp/build/include
- Double check cuDecomp building is fine (must be compiled using HPC-SDK)
- Single folder: contains the single GPU version of the code (see MHIT36 repository for further details), no MPI required.
- Multi folder: multi GPU version of the code (work in progress). Use local.sh or leo.sh to compile and run the code (see porting status above); the multi GPU version relies on cuDecomp for Pencil Transposition and halo exchanges.
- Autotuning of the multi-GPY version: leave pr=0 and pc=0, cuDecomp will perform an autotuning at the start finding the best decomposition (the only input is the total number of tasks). Everything should be automatic in the code (as it is obtained from cuDecomp variable)