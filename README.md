# Pervaporation

In the file "membrane_module" there is a file named "run_this_file", you should run this fil; it is the main file.

The spreadsheets have values, coefficients, data points, which the program uses to model water and ethanol in a pervaporation unit. You can change these values if you wish.

The model is based on a paper by Holger Thiess, Axel Schmidt, and Jochen Strube (DOI: 10.3390, also found here: https://doi.org/10.3390/membranes8010004). It accounts for concentration polarisation at the membrane surface, the effects of water on permeance, and general diffusion through non-porous membranes. Be careful with the study though, as they model concentration gradients across the width of the membrane / perpendicular to the membrane surface. This program models the change in weight fractions as a stream passes through a pervaporation unit, assuming two rectangular membranes on each side of the flow.

The utility files contain class objects for calculating thermodynamic and stream properties. I did this to streamline the ODEs, otherwise it'd become a tangled mess of equations -- no thanks.

It is very basic, and it served its purpose well.

One caveat, the visocsity of the stream wasn't calculated and was left blank. You can add the Erlying method if you wish to do so.

For further information, send me an email at benj.alexander.cox@proton.me
