% This system is of advanced type
% Not solvable with our solver

clear all; close all; clc

E=@(t)[t 0; 0 0];
A=@(t)[0 0; t 0];
B=@(t)[0 1; 1 1];

%phi = @(t) [sin(t);cos(t)];
%dphi = @(t)[cos(t);-sin(t)];
phi = @(t) [exp(t/10);t];
dphi = @(t)[0.1 * exp(t/10);1];
    
tau = @(t) exp(-t);

f=@(t)E(t)*dphi(t)-A(t)*phi(t)-B(t)*phi(t-tau(t));