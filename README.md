
# Moleculer Connectivity Analysis

## Table of contents
* [Aim](#Aim)
* [Method](#Method)
* [Technologies](#technologies)
* [Setup](#setup)
* [Status](#Status)

## Aim 

the aim of this work is to develop a simple, general and robust automated algorithm for
determining connectivities and topologies of a wide variety of molecular systems.

## Method

Automatically detecting molecular connectivities and topologies, without requiring user input or defining system-specific sets of rules, is a challenging algorithmic problem. The novel strategy developed in this work relies on two overlapping, or near equivalent, assumptions; one implicit and one explicit. The explicit assumption is that the closer two atoms are the more likely they are to be involved in a bonding interaction. The second implicit assumption is that the atomic positions are chemically sensible and exist in a stable "chemically reasonable" configuration where bonds are not actually in flux.

This work also relies on the novel insight that it should be simpler to create more bonds than necessary and remove those that are not physically meaningful, rather than to pre-define a rigid set of rules that only allows physically meaningful bonds to be created in the first place. Similarly, it should be easier to find rings by removing non-ring components than pre-defining rules to identify them.

In brief, there are four key stages to determining molecular connectivities and topologies:
1. Read in and store atomic coordinates and identities
2. Process atomic coordinate data to find all physically meaningful bonds
3. Process bond connectivity information to identify connected ring systems
4. Assign and store formal bond types and charges.


## Technologies 

This project is written in Python Version 3.5.4

libraries included:

* force_field_code_final
* ringalgothem_final_2
* MD_main_final_2  
* expected_ring_number 
* time
* numpy
* math
* re
* pdb_formating_data 
* ringalgothem_final_2 
* mpl_toolkits.mplot3d
* matplotlib.pyplot
* collections
* periodic_table 

## Setup

To run this project, download the files locally into the same directory:

run organic_pdb_processing_final.py 

## Status # rewrite

The project was complete in 2018 achiving both of the major goal. Which where
to construct an automated method for constructing a set of 3D connectivities
for a molecule with a reasonable set of Cartesian coordinates and a automated
method for finding and identifing ring systems within a molecule for which 
there is a complete and constant set of 3D connectivities.

A in-depth descriptions of the above methods can be found in my masters thesis 
which can be found here: https://ir.canterbury.ac.nz/handle/10092/16823
