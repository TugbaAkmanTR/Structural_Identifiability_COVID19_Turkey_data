%% Full model 

close all; clear all; clc;

set(0, 'defaultaxesfontsize',16)

% initialize the model
Setup_StageI

%% Optimization
arSimu(true,true,true);

n_fit = 500; % number of starts

if n_fit > 1
    arFitLHS(n_fit); % multistart optimization
else
    arFit;
end

% Print optimization results
arPrint;     % display parameter values
arSimu(true,true,true);

arPlot

%arSave('Model')

%Profile likelihood
arPLEInit
ple


%% To reproduce the results in the paper
% arLoad('20241206T152503_Model')
% arPlot
%arPLEInit
%ple
