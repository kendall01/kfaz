function [ MFP ] = findMFP( Ma,k )
MFP = Ma*sqrt(k)*( 1 + (Ma^2)*( (k-1)/2) )^(-(k+1)/(2*(k-1)));
end

