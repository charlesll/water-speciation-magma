# water-speciation-magma

Model for prediction of water speciation in magmatic liquids

# Citation

This software was part of the publication:

Moretti, R., Le Losq, C., and Neuville, D.R. (2014) The amphoteric behavior of water in silicate melts from the point of view of their ionic-polymeric constitution. Chemical Geology, 367, 23–33. https://doi.org/10.1016/j.chemgeo.2013.12.012

This should be cited in case of use.

# Licence

MIT licence, see Licence file.

# Contributors

Roberto Moretti (IPGP), moretti@ipgp.fr

Charles Le Losq (IPGP), lelosq@ipgp.fr

# Dependencies

A working fortran compiler. We suggest using gfortran, tested on Mac and Linux. It works well with this software!

# How to use

Download the repository, and use the provided example input file. It first requires compilation of the FORTRAN source, then running the compilated software.

## Compilation

To create the program, with gfortran, in the terminal on Linux or MacOS:

`$ gfortran waspecCLL.f90 -o waspec.o`

## Running the software

The software takes an input file, waspec.in, which contains the compositions of interest.

It returns an output file, waspec.out

Run in the terminal, after compilation, run the command:

`$ ./waspec.o`

## How to use the input file waspec.in

1) first line is an integer (number of compositions) followed by another integer, the KYMPA variable (1 or 0);
this is explained later on

4) write the composition line in the form of:

T(°C) P(bar) SiO2 TiO2 Al2O3 Fe2O3 Cr2O3 FeO Mno MgO Cao Na2O K2O P2O5 H2Otot H2Omol

NOTES:

- Oxides are given in wt%.
- H2Omol is the molecular water in wt% (otherwise put a number, such as 0.1, in any case lower than H2Otot)
- T(°C) is the temperature of the liquid (greater than or equal to Tg)

## About waspec.out:

1) If you have set to 0 the second integer, right after the number of composition, waspec considers input H2Omol as an experimental constraint. If you have set it to 1, waspec iterates until convergence. Therefore, results can be different between the "0" and "1" KYMPA options.

2) The output lists many parameters:

- H2Ototwt%: total water in wt%, given in input.
- nH2Otot: total moles of water per 100g of melt (inout value).
- nhtot: total moles of hydrogen per 100g of melt (input value).
- nH2Omolwt%EXP: molecular water in wt% (input value).
- nH2OmolEXP: moles of molecular water per 100g of melt (input value).
- 1/T: the reciprocal of temperature (in Kelvin).
- nfreeOH: calculated moles of OH- per 100g of melt.
- nfreeH: calculated moles of H+ per 100g of melt.
- nOH-IR: the sum of calculated moles of H+ and OH- per 100g of melt.
- nH2OmolEXP: moles of molecular water per 100g of melt (input value) (same as above).
- nH2Omolcalc: calculated moles of molecular water per 100g of melt.
- nHtot_after: calculated moles of total hydrogen per 100g of melt (CAN BE DIFFERENT THAN nhtot, IF KYMPA = 0)
- DISPERSION: the difference between nhtot_after and nhtot
- XfreeOHcalc: mol fraction of calculated free OH
- XfreeHcalc: mol fraction of calculated H
- XH2Omolcalc: mol fraction of calculated H2Omol
- YfreeOHcalc: relative proportion of calculated free OH
- YfreeHcalc: relative proportion of calculated free H
- YH2Omolcalc: relative proportion of calculated H2Omol

