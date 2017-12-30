# SHG FROG software for measuring n2 with bootstrap FROG error

This program when finished is supposed to get measured SHG FROG spectrogram as input and give n2 and n2 error as output.

Detailed functioning:
- get Femtosoft FROG compatible SHG FROG spectrogram matrix as input
- retrieval of the complex pulse through (principal component) generalized projections 
- computing FROG error with bootstrap method
- fitting experimental data through Total Least Squares method including all possible sources of error
- get n2 and its error as output
