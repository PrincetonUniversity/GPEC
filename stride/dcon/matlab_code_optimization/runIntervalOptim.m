%This code is run to optimize the estimation of runtime given
%certain distances from singular points.

%The output y houses the desired weightings, comparing two different
%models (with 200 and 300 intervals, e.g.)

clear all;
close all;
load('intervalOptimData.mat');
x = getIntervalOptim(psi12, s1to3, t12);
load('intervalOptimData2.mat');
x2 = getIntervalOptim(psi12, s1to3, t12);

y = [x x2];  %Look at y in the Workspace!