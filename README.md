# FFT_bootstrap

This is the code that I generally use in my trend emergence stuff.
This will be getting updated to make it more generally useful, but 
for now it should work out of the tin.

Currently, it should work for either RAPID or NEMO data. It's been
updated to make it significantly faster, but as a result is using a 
pared down polynomial fitting function. This might cause issues if 
the timeseries being used are dodgy.
