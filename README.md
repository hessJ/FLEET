# FLEET
Functional LD-interval EnrichmEnt Test (active development, i.e., easy to break)

This tool is used for biological interpretation of GWAS signals for polygenic phenotypes. 

### Please visit the Wiki page to set up FLEET on your workstation
https://github.com/hessJ/FLEET/wiki

#### Download and install

     git clone https://github.com/hessJ/FLEET

#### Check for updates

     git pull

A description of the method is provided in the powerpoint slides (fleet_slides.ppt).

## Run 

    bash install.sh

## Display help commands/options:

    fleet.R --help 
    	or
    fleet.R -h
  
  Will print something like this: 
  
  ```  
 
==========================================================
*
* Functional LD-interval EnrichmEnt Test (FLEET)
*
* Jonathan L. Hess, PhD and Stephen J. Glatt, PhD (c) 2017
*
* SUNY Upstate Medical University, PsychGENe Lab
*
* Contact: hessjo@upstate.edu
*
* https://github.com/hessJ/FLEET
*
* GNU GENERAL PUBLIC LICENSE v3
===========================================================

Start time: 2017-10-04 11:25:16

Location of fleet: /Users/jonathanhess/Documents/FLEET/fleet.R
Usage: /Users/jonathanhess/Documents/FLEET/fleet.R [options]


Options:
	-G CHARACTER, --gwas=CHARACTER
		Path to GWAS summary statistics. Column headers are required. Allowed delim = sep, tab, or comma

	-O CHARACTER, --out=CHARACTER
		output file name [default = fleetOut]

	-R DOUBLE, --r2=DOUBLE
		R-squared threshold for linkage disequilibrium calculations [default = 0.6]

	-W INTEGER, --ld-window=INTEGER
		Size of window (kilobases) for calculating linkage disequilibrium [default = 1000]

	-S CHARACTER, --snp-field=CHARACTER
		SNP column header in GWAS file

	-P CHARACTER, --pcol=CHARACTER
		P-value column header in GWAS file

	--robust=DOUBLE
		Computing robust standard errors using White method (via vcovHC function in sandwich pkg) [default = TRUE]

	-N INTEGER, --nPerms=INTEGER
		Number of permutations to perform [default = 1000]

	-T INTEGER, --threads=INTEGER
		Number of cores for parallel operations [default = 1]

	-L CHARACTER, --label-annotations=CHARACTER
		Path to annotation table [default = /Users/jonathanhess/Documents/FLEET/annotations/annotation.txt]

	-D CHARACTER, --rd-annots=CHARACTER
		Path to .Rdata annotations [default = /Users/jonathanhess/Documents/FLEET/annotations/]

	-F DOUBLE, --annot-cnt=DOUBLE
		Minimum annotation count observed across LD-clumps [default = 10]

	-M LOGICAL, --fleet-prune-ref=LOGICAL
		Initiate LD-pruning step of 1KG reference data. Only needs to be run once. [default = TRUE]

	-A LOGICAL, --fleet-annotate=LOGICAL
		Annotate LD-clumps with bedtools [default = TRUE]

	-E LOGICAL, --fleet-enrichment=LOGICAL
		Perform enrichment analysis with weighted linear models [default = TRUE]

	--fleet-permutation=LOGICAL
		Perform enrichment analysis with permutation (randomizing annotations) [default = TRUE]

	--fast-permutation=LOGICAL
		Simple permutation analysis [default = FALSE]

	--robust-permutation=LOGICAL
		Permutation analysis that will sample variants from the MAF bin of target SNPs [default = FALSE]

	--speed=CHARACTER
		Change behavior of linear models (fast mode: SET becomes response variable, slow mode: Z-score becomes response variable) [default = fast]

	--pthres=CHARACTER
		Table with P-value threshold(s) for SNP bins [default P-values < 5e-08, 1e-07, and 1e-06]

	--plots=LOGICAL
		Turning this on will produce multiple plots to display summary statistics [default = FALSE]

	-h, --help
		Show this help message and exit

```




## GNU GENERAL PUBLIC LICENSE v3

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

========================================================================

Made on a MacBook Pro Retina, 2.5 GHz Intel Core i7, 16 GB 1600 MHz DDR3
