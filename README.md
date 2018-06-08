# Polygon400 grid generator for sCRACM

This script generates a 12x24 grid of stimuli (228 patterns) to use as input 
for the Mightex Polygon400. 

Each pattern contains one stimulus spot, in the form of an array of 0s with 1s at the spot location.
Different spot locations are separated by one line of 0s. 
The size of the array thus determines the grid resolution and the relative sizes of the
light spots and spaces.

After generating an ordered array of stim patterns, the order of the patterns 
is randomly shuffled and the order with the greatest separation between 
subsequent patterns is chosen. This array is the exported to a .txt file that
can be read by the Polygon400.
