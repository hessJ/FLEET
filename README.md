# FLEET
Functional LD-clump EnrichmEnt Test 

This tool is used for biological interpretation of GWAS signals for polygenic traits. 

### Please visit the Wiki page to set up FLEET on your workstation
https://github.com/hessJ/FLEET/wiki


A description of the method is provided in the powerpoint slides (fleet_slides.ppt).

## FLEET is a modular script

1. `LD pruning` in Module A
2. `GWAS clumping` in Module B
3. `Annotating hg19 variants` in Module C
4. `Enrichment tests` in Module D
5. `Plotting and exporting html report` in Module E

### Tip for first-time use:
1. Run `Module A` once to get LD-pruned 1000G data
2. Skip `Module A` in future runs if satisfied with the LD-pruning parameters
