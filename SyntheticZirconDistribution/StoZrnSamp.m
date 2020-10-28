
function [out] = StoZrnSamp(sim,mod,nIt,nZr)

% StoZrnSamp: randomly draws a selected number of relative zircon
%              crystallization ages from a synthetic cumulative zircon
%              crystallizatioon distribution. The drawing procedure can 
%              be iterated an user defined number of times.
%
% Arguments: (input)
%   sim - Integer (1 to 24). The user can select between 24 rhyolite-MELTS 
%         based synthetic zircon crystallization distributions. The variable 
%         values (1-24) are related to different simulation scenarios, which  
%         details are contained in Tab. DR4 of the supplementary material of 
%         Tavazzani et al., 2020
%
%   mod - Inyteger (0 or 1). The user can select outputs from either 
%         plutonic (mod=0) or volcanic (mod=1) crystallization simulations.
%          
%   nIt - Integer (1 to n). Number of Iterations of the random sampling
%         process (to ensure a robust result nIterations should be >= 1000).
%
%   nZr - Integer (1 to n). Number of stochastically sampled zircon ages 
%         from choosen synthetic zircon crystallization probability 
%         distribution.For comparison of synthetic zircon age spectra with
%         empirical spectra from a natural sample, nZir should be adjusted
%         based on the number of zircons dates obtained on the natural
%         sample of interest.
%  
%
% Author: Lorenzo Tavazzani
% e-mail address: ltavazzani@smu.edu
% Release: 1.0
% Release date: 10/23/20

% Load synthetic zircon crystallization distributions
t=xlsread('dZr_Matlab_set2.xlsx',1); %dZrt
p=xlsread('dZr_Matlab_set2.xlsx',2); %dZr
tc=xlsread('dZr_Matlab_set2_cutoff.xlsx',1); %dZrt_cutoff at 740°C after last recharge (eruption)
pc=xlsread('dZr_Matlab_set2_cutoff.xlsx',2); %dZr_cutoff at 740°C after last recharge (eruption)%Select P (0 or 1) at the beginning of routine to use cutoff temperature

% Set mode of crystallization (Plutonic vs. Volcanic)
if mod==1
    t=tc; p=pc;
else
    t=t; p=p;
end

% Select simulation of interest and eliminates NaN values from matrices of 
% zircon crystallization probabilities 
p=p(:,sim);
t=t(:,sim)';
p(isnan(p)) = [];
t(isnan(t)) = [];

% Smoothing and normalization of zircon crsytallization probabilities to 1
p=smooth(p)';
p=p./sum(p);  

% Randomly generates a population of synthetic zircon crystallization ages
% using the Matlab Built-in function randsrc()
out = randsrc(nIt,nZr,[t;p]);

% Order simulated zircon dates "youngest" to "oldest"
a=size(out); a=a(:,1);
for ii=1:a
out(ii,:)=sort(out(ii,:),'ascend');
end

% Absolute zircon crystallization times are scaled between onset of zircon
% saturation and crystallization of youngest zircon 
out=out';
for i=1:size(out,2)
    out(:,i) = out(:,i) - min(out(:,i));
    out(:,i) = out(:,i)./max(out(:,i));
end

end