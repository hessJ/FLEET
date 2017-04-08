# FLEET
Functional LD-clump EnrichmEnt Test (active development)

This tool is used for biological interpretation of GWAS signals for polygenic traits. 

### Please visit the Wiki page to set up FLEET on your workstation
https://github.com/hessJ/FLEET/wiki


A description of the method is provided in the powerpoint slides (fleet_slides.ppt).

## FLEET is a modular script
The first two modules are wrappers for Plink. The last three modules perform the main operations of FLEET and produce output.

1. `LD pruning` in Module A (Spares genome-wide significant markers)
2. `GWAS clumping` in Module B
3. `Annotating hg19 variants` in Module C
4. `Enrichment tests` in Module D
5. `Plotting and exporting html report` in Module E

This script is under active development. Stay tuned for updates.

### Tip for first-time use:
1. Run `Module A` once to store LD-pruned 1000G data in `~/PATH_TO_FLEET/pruned/` directory
2. Skip `Module A` in future runs if satisfied with the LD-pruning parameters
