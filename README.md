# FFT_bootstrap

This is the code that I generally use in my trend emergence stuff.
This will be getting updated to make it more generally useful, but 
for now it should work out of the tin.

Currently, it should work for either RAPID or NEMO data. It's been
updated to make it significantly faster, but as a result is using a 
pared down polynomial fitting function. This might cause issues if 
the timeseries being used are dodgy.

Files titles FFT_*_CI_legacySlow.m are using the older, safer, slower 
code. If using a new dataset with this code, I would recommend
checking consistency between the two before doing anything intensive.
