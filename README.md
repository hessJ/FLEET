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

## Display help commands/options:

    Rscript FLEET.R -h
  
  Will print this: 
  
  ```  
  ==========================================================
  *
  * Functional LD-clump EnrichmEnt Test (FLEET)
  *
  * Jonathan L. Hess, PhD and Stephen J. Glatt, PhD (c) 2017
  *
  * SUNY Upstate Medical University, PsychGENe Lab
  *
  * Contact: hessjo@upstate.edu
  *
  * GNU GENERAL PUBLIC LICENSE v3
  ===========================================================
Usage: FLEET.R [options]


Options:
	-g CHARACTER, --gwas=CHARACTER
		Path to GWAS summary statistics

	-o CHARACTER, --out=CHARACTER
		output file name [default = out.txt]

	-r2 DOUBLE, --r2=DOUBLE
		R-squared threshold for linkage disequilibrium calculations [default = 0.2]

	-ldw INTEGER, --ld-window=INTEGER
		Size of window (kilobases) for calculating linkage disequilibrium [default = 1000]

	-s CHARACTER, --snp-field=CHARACTER
		SNP column header in GWAS file

	-c CHARACTER, --clump-field=CHARACTER
		P-value column header in GWAS file

	-p DOUBLE, --permStartThreshold=DOUBLE
		Minimum P-value from enrichment test to initiate permutation analysis [default = 0.005]

	-n INTEGER, --nPerms=INTEGER
		Number of permutations to perform [default = 1000]

	-t INTEGER, --threads=INTEGER
		Number of cores for parallelization [default = 1]

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

Made on a 
	MacBook Pro Retina
	2.5 GHz Intel Core i7
	16 GB 1600 MHz DDR3
