# SHG FROG PCGPA software for measuring nonlinear refractive index

This program will take measured SHG FROG spectrogram and some experimental parameters as input and give n2 with error as output for particular measured material.

FROG algorithm with bootstrap implemented is working correctly in MATLAB.
After everything works properly in matlab it will be rewritten in Python in the future.

Main part of this software (FROG algorithm) is based on:

- Adam Wyatt's program, available here: https://www.mathworks.com/matlabcentral/fileexchange/16235-frequency-resolved-optical-gating--frog-
- Steven Byrnes' program (which is direct update of the first one), available here: https://www.mathworks.com/matlabcentral/fileexchange/34986-frequency-resolved-optical-gating--frog-
- Kenneth DeLong's program, available here: https://github.com/kenwdelong/frog
- many publications by Rick Trebino, Daniel Kane, Kenneth DeLong, Ziyang Wang et al.
