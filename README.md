# FLEET
Functional LD-clump EnrichmEnt Test (active development, i.e., easy to break)

This tool is used for biological interpretation of GWAS signals for polygenic phenotypes. 

### Please visit the Wiki page to set up FLEET on your workstation
https://github.com/hessJ/FLEET/wiki

#### Download and install

     git clone https://github.com/hessJ/FLEET

#### Check for updates

     git pull

A description of the method is provided in the powerpoint slides (fleet_slides.ppt).

## FLEET is a modular script
The first two modules are wrappers for Plink. The last three modules perform the main operations of FLEET and produce output.

1. `LD pruning` in Module A (Spares genome-wide significant markers)
2. `GWAS clumping` in Module B
3. `Annotating hg19 variants` in Module C
4. `Enrichment tests` in Module D
5. `Plotting and exporting html report` in Module E

### Tip for first-time use:
1. Run `Module A` once to store LD-pruned 1000G data in `~/PATH_TO_FLEET/pruned/` directory
2. Skip `Module A` in future runs if satisfied with the LD-pruning parameters


#### GNU GENERAL PUBLIC LICENSE v3

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
