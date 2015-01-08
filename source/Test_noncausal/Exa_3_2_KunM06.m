% This example is the modified version of Example 3.2 in the blue book of
% Volker about DAEs.

clear all; close all; clc

tau = @(t)1;

E=@(t)[0 0; 1 -t];
A=@(t)[0 0; 0 0];
B=@(t)[1 -t-tau(t);0 0];

%phi = @(t) [sin(t);cos(t)];
%dphi = @(t)[cos(t);-sin(t)];
phi = @(t) [exp(t/10);t];
dphi = @(t)[0.1 * exp(t/10);1];
    

f=@(t)E(t)*dphi(t)-A(t)*phi(t)-B(t)*phi(t-tau(t));

xe = phi;
