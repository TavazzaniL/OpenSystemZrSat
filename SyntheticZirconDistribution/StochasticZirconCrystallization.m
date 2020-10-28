%% Randomly generates zircon age distributions based on MCS outputs of 
% relative zircon growth under different crystallization scenarios
%  
% Author: Lorenzo Tavazzani
% e-mail address: ltavazzani@smu.edu
% Release: 1.0
% Release date: 10/23/20

clear all
close all

                          %%%%%% [User Input] %%%%%% 
% Select path to folder that contain the scripts and .xls files
cd /Users/Utente/Desktop/SyntheticZirconDistribution

                          %%%%%% [User Input] %%%%%% 
%%%%%% Parameters used for stochastic sampling of  [User Input] %%%%%%
% Stochastic Zircon sampling function, user adjustable parameters 
% [out] = StoZrnSamp(mod,sim,nZr,nSim)
% (see StoZrnSamp.m for function usage and parameters description) 
[out] = StoZrnSamp(1,0,100,100);

                          %%%%%% [User Input] %%%%%% 
%%%%%% Parameters used for the graphical display of data [User Input] %%%%%%
% [nZr] Set the number of randomly selected zircons per iteration, it needs
% to be adjusted based on the nZr parameter used to run the function StoZircSamp
nZr = 100; 
% [alpha] Set the transparency of the multiple outputs from the stochastic 
% sampling process displayed as KDEs. Need to be adjusted based on the
% nSim (alpha = 0.2 for 100 < nSim < 1000) (alpha = 0.08 for nSim > 1000)
alpha = 0.2; 

                          %%%%%% [User Input] %%%%%% 
%%%%%% Load and plot empirical zircon age distribution for visual comparison
% Select path to folder that contain the .csv files of empirical zircon 
% ages distribution and the Matlab script IDTIMS_DataLoad_KDE 
cd /Users/Utente/Desktop/LiteratureDatasets   
IDTIMS_DataLoad_KDE
% Select up to four ID-TIMS dataset to compare with modeling simulation results
% References stored in IDTIMSreferences.txt  
% (**Intrusive bodies**: Valle_Mosso_Pluton / Valle_Mosso_Pluton_NoAntecrysts 
%      Bergell_pluton / Capanne_pluton / Mount_Givens_pluton / Tenpeak_pluton  
%      Mount_Stuart_pluton / JohnMuir_IntrusiveSuite / Tuolomne_IntruiveSuite
%      LagodellaVacca_IntrusiveComplex / ValFredda_IntrusiveComplex 
%      LaGloria_pluton / Golden_Horn_Batholith / Turkey_Creek_Intrusives
%      Latir_Intrusives)
% (**Volcanic eruptions**: Sesia_Caldera / Sesia_Caldera_NoAntecrysts  
%      Kneeling_Nun_Tuff / Huckleberry_Ridge_Tuff / Lava_Creek_Tuff 
%      Fish_Canyon_Tuff / Turkey_Creek_Caldera / Latir_Volcanics)
% Insert name of the igneous unit from the list above (check spelling)
EmpDist_1=Bergell_pluton; 
EmpDist_2=Capanne_pluton;
EmpDist_3=JohnMuir_IntrusiveSuite;
EmpDist_4=Golden_Horn_Batholith;
% Insert igneous unit name for legend visualization (spelling does not matter)                       
EmpDist_1_str='Bergell pluton'; 
EmpDist_2_str='Capanne pluton'; 
EmpDist_3_str='John Muir Intrusive Suite';
EmpDist_4_str='Golden Horn Batholith';



                   %%%%%% [User editing not necessary] %%%%%%   
%%%%%% Graphical output   
% Kernel Density Estimate of averaged zircon age spectra, scaled between 
% initiation and termination of zircon crystallization
[out_f,out_x,out_bw] = ksdensity(mean(out, 2));
outdistAvg = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,100));
x = linspace(0,1,100); outdistAvg = outdistAvg./trapz(x,outdistAvg);

% Plot average Kernel Density Estimate
figure(1)
plot(x,outdistAvg,'Color','k','LineStyle','-','Linewidth',2)
hold on; 

% Plot zircon age empirical datasets KDEs
hold on; 
c = winter(4);
plot(x,EmpDist_1,'Color',c(1,:),'LineStyle','--','Linewidth',2)
plot(x,EmpDist_2,'Color',c(2,:),'LineStyle','-','Linewidth',2)
plot(x,EmpDist_3,'Color',c(3,:),'LineStyle','-.','Linewidth',2)
plot(x,EmpDist_4,'Color',c(4,:),'LineStyle',':','Linewidth',2)
legend('Synthetic Zrn (avg.)',EmpDist_1_str,EmpDist_2_str,EmpDist_3_str,EmpDist_4_str,'AutoUpdate','off');

% Kernel Density Estimates of all zircon age spectra, scaled between 
% initiation and termination of zircon crystallization
for ii=1:nZr
   [out_f,out_x,out_bw] = ksdensity(out(:,ii));
   outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,100));
   x = linspace(0,1,100); outdist = outdist./trapz(x,outdist); 
   % Plot all Kernel Density Estimates
   plot1=plot(x,outdist,'Color','k','LineStyle','-','Linewidth',0.15);
   plot1.Color(4) = alpha;
end

% Plot details
ylim([0 3]);
set(gca,'xdir','reverse')
xlabel("Age (scaled)")
ylabel("Density")
box on
