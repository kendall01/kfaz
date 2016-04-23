%pemStoich.m
%4-22-16 Created Jon Renslo
function [ mixVec ] = pemStoich( lambda, mH2)
%PEMSTOICH 
% does stoichiometry calculation for a PEM H2 and air fuel cell
% mH2 in kg
% mixVec = [ water vapor, liquid water, o2, n2] in moles of product formed
MMH2 = 2.02; %g/mol
gToKg = 1000;

Nh2 = mH2/MMH2*gToKg;

%TODO based on exit conditions, must calculate liquid/gas water phase for products

mixVec = [0];


end