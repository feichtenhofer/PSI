%Reference:
%   
%   C Feichtenhofer, H Fassold, P Schallauer
%   "A perceptual image sharpness metric based on local edge gradient
%   analysis", IEEE Signal Processing Letters, 20 (4), 379-382, 2013
%   
%
%   Written by Christoph Feichtenhofer (cfeichtenhofer AT gmail.com)
%   feichtenhofer.github.io                           
%

load mandrill; I=X;

psi =  PSI(I)
figure;im(I)