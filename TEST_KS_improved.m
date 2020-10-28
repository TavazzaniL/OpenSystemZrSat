%% Randomly generated zircon samples based on MCS outputs of relative zircon 
% growth under different crystallization scenarios 
% 1/7/2020 Lorenzo Tavazzani

clear all
close all
clc

% Data from literature used in bootstrapped crystallization dist. plots
cd /Users/Utente/Desktop/zurich_ID-TIMS/Keller_papers_zircSat/BayeZirChron.c-master/examples/BootstrappedCrystallizationDistributionExamples/SCRIPTS
IDTIMS_DataLoad

%% Randomly generated samples selected zircon crystallization probability denisty function
% (sample PST214 composition)
close all
clear all
clc

% Parameters
c = lines(7);
nSimulations=100; 
alpha=0.2; % Adjust transparency based on number of simulations (nSimulations)
Alpha=0.2; % Significance level of KS test
nZircons=30; % Number of stochastically sampled zircons

P=0; %Volcanic (P=1) or plutonic (P=0) crystallization simulations
S=14; % ID of simulation on Tab. DR4 of supplementary material

%Loading simulations datasets
cd /Users/Utente/Desktop/zurich_ID-TIMS/4_Modeling/Zr_sat_runs/1_MCS_runs/1_dZr_Results
t=xlsread('dZr_Matlab_set2.xlsx',1); %dZrt
p=xlsread('dZr_Matlab_set2.xlsx',2); %dZr
tc=xlsread('dZr_Matlab_set2_cutoff.xlsx',1); %dZrt_cutoff at 740°C after last recharge (eruption)
pc=xlsread('dZr_Matlab_set2_cutoff.xlsx',2); %dZr_cutoff at 740°C after last recharge (eruption)%Select P (0 or 1) at the beginning of routine to use cutoff temperature

if P==1
    t=tc; p=pc;
else
    t=t; p=p;
end

p=p(:,S);
t=t(:,S)';
p(isnan(p)) = [];
t(isnan(t)) = [];
p=smooth(p)'; % smoothing filter on the data, can be avoided
p=p./sum(p);  % normalization of growth rates to p=1

% randomly generated population
out = randsrc(nSimulations,nZircons,[t;p]);

% order "youngest" to "oldest"
a=size(out); a=a(:,1);
for ii=1:a
out(ii,:)=sort(out(ii,:),'ascend');
end

figure(2)
% histogram of randomly generated samples
out=out';

b=size(t,2); 
hold on
for i=1:size(out,2)
   %histogram(out(:,i),linspace(min(t),max(t),max(t)-min(t)+1)); hold on
end

xlim([min(t)-1 max(t)+1]); ylim([0 nZircons/2])
set(gca,'xdir','reverse'); box on
xlabel("t interval"); ylabel("N"); 

% Plot (rank order) data
plot(out(:),'.','MarkerSize',10)
ylim([min(t) max(t)])
xlabel("N"); ylabel("Age (No unit)")

% Plot (rank order) scaled data
for i=1:size(out,2)
    out(:,i) = out(:,i) - min(out(:,i));
    out(:,i) = out(:,i)./max(out(:,i));
end
plot(out(:),'.','MarkerSize',10)
xlabel("N"); ylabel("Age, scaled")

%%%%%%%%%
%Kolmogorov-Smirnov test setup for 1 random simulation

subplot(2,2,1)
cd /Users/Utente/Desktop/zurich_ID-TIMS/Keller_papers_zircSat/BayeZirChron.c-master/examples/BootstrappedCrystallizationDistributionExamples/UPB_LIT
bergell=csvread('BergellAll.csv'); bergell(bergell==0)=NaN;
for i=1:size(bergell,2)
    bergell(:,i) = bergell(:,i) - min(bergell(:,i));
    bergell(:,i) = bergell(:,i)./max(bergell(:,i));
end
nn=round(rand.*100); nb=2; ns=4;
a=out(:,nn); b=bergell(:,nb);
[h,p] = kstest2(a,b,'Alpha',Alpha) %Kolmogorov-Smirnov test

%%%%%%%%%
SesiaG=csvread('SesiaG.csv');
SesiaG(SesiaG==0)=NaN;
for i=1:size(SesiaG,2)
    SesiaG(:,i) = SesiaG(:,i) - min(SesiaG(:,i));
    SesiaG(:,i) = SesiaG(:,i)./max(SesiaG(:,i));
end
e=SesiaG(:,ns);
[h1,p1] = kstest2(a,e,'Alpha',Alpha) %Kolmogorov-Smirnov test

hold on; ecdf(a); ecdf(b); ecdf(e);
set(gca,'xdir','reverse'); set(gca,'ydir','reverse'), box on
legend('Synthetic Zrn','Bergell','SesiaG','Location','southeast');

%%%%%% Bootstrapped KDEs
npoints = 100;
%[out_f,out_x,out_bw] = ksdensity(mean(out, 2));
[outKS_f,outKS_x,outKS_bw] = ksdensity(out(:,nn));
[outBKS_f,outBKS_x,outBKS_bw] = ksdensity(bergell(:,nb));
[outSGKS_f,outSGKS_x,outSGKS_bw] = ksdensity(SesiaG(:,ns));
%outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,npoints));
outdistKS = interp1(outKS_x,outKS_f,linspace(0-0.05,1+outKS_bw,npoints));
outdistBKS = interp1(outBKS_x,outBKS_f,linspace(0-0.05,1+outBKS_bw,npoints));
outdistSGKS = interp1(outSGKS_x,outSGKS_f,linspace(0-0.05,1+outSGKS_bw,npoints));
x = linspace(0,1,npoints);
%outdist = outdist./trapz(x,outdist);
outdistKS = outdistKS./trapz(x,outdistKS);
outdistBKS = outdistBKS./trapz(x,outdistBKS);
outdistSGKS = outdistSGKS./trapz(x,outdistSGKS);

%%%%% Plotting KDEs
subplot(2,2,2)
hold on; 
%plot(x,outdist,'Color','k','LineStyle','-','Linewidth',2)
plot(x,outdistKS,'Color',c(1,:),'LineStyle','-','Linewidth',1)
plot(x,outdistBKS,'Color',c(2,:),'LineStyle','-','Linewidth',1)
plot(x,outdistSGKS,'Color',[0.9290, 0.6940, 0.1250],'LineStyle','-','Linewidth',1)

legend('Synthetic Zrn',['p= ', num2str(round(p,3))],['p= ', num2str(round(p1,3))],'Location','northwest');
ylim([0 3]);
set(gca,'xdir','reverse')
xlabel("Age (scaled)")
ylabel("Density")
box on

for ii=1:nSimulations
a=out(:,ii); 
[h,p] = kstest2(a,b,'Alpha',Alpha); %Kolmogorov-Smirnov test
hresult(ii) = h;
presult(ii) = p; presult=sort(presult,'descend');
[h1,p1] = kstest2(a,e,'Alpha',Alpha); %Kolmogorov-Smirnov test
h1result(ii) = h1;
p1result(ii) = p1; p1result=sort(p1result,'descend');
end

ptest=size(find(hresult==0),2)
p1test=size(find(h1result==0),2)

%
subplot(2,2,3)
hold on; 
plot(x,outdistSGKS,'Color','k','LineStyle','-','Linewidth',1)
outpass=out(:,find(h1result==0)); % percent passed KS compared to VMP
%plot(x,outdistBKS,'Color','k','LineStyle','-','Linewidth',1)
%outpass=out(:,find(hresult==0)); % percent passed KS compared to Bergell
for ii=1:size(outpass,2)
   [outpass_f,outpass_x,outpass_bw] = ksdensity(outpass(:,ii));
   outpassdist = interp1(outpass_x,outpass_f,linspace(0-0.05,1+outpass_bw,npoints));
   x = linspace(0,1,npoints);
   outpassdist = outpassdist./trapz(x,outpassdist); 
   plot1=plot(x,outpassdist,'Color',c(1,:),'LineStyle','-','Linewidth',1);
   plot1.Color(4) = alpha;
end

x1= 0.98;
y1= 2.8;
text(x1,y1,['KS (p>',num2str(Alpha),') = ', num2str(p1test),'%'],'FontSize',10,'Color',c(1,:));
%text(x1,y1,['KS (p>',num2str(Alpha),') = ', num2str(ptest),'%'],'FontSize',10,'Color',c(1,:));

ylim([0 3]);
set(gca,'xdir','reverse')
xlabel("Age (scaled)")
ylabel("Density")
box on

subplot(2,2,4)
hold on;
plot(x,outdistSGKS,'Color','k','LineStyle','-','Linewidth',1)
outfail=out(:,find(h1result==1)); % percent passed KS compared to VMP
%plot(x,outdistBKS,'Color','k','LineStyle','-','Linewidth',1)
%outfail=out(:,find(hresult==1)); % percent passed KS compared to Bergell
for ii=1:size(outfail,2)
   [outfail_f,outfail_x,outfail_bw] = ksdensity(outfail(:,ii));
   outfaildist = interp1(outfail_x,outfail_f,linspace(0-0.05,1+outfail_bw,npoints));
   x = linspace(0,1,npoints);
   outfaildist = outfaildist./trapz(x,outfaildist); 
   plot1=plot(x,outfaildist,'Color','r','LineStyle','-','Linewidth',1);
   plot1.Color(4) = alpha;
end

x2= 0.98;
y2= 2.8;
text(x2,y2,['KS (p<',num2str(Alpha),') = ',num2str(100-p1test),'%'],'FontSize',10,'Color','r');
%text(x2,y2,['KS (p<',num2str(Alpha),') = ',num2str(100-ptest),'%'],'FontSize',10,'Color','r');

ylim([0 3]);
set(gca,'xdir','reverse')
xlabel("Age (scaled)")
ylabel("Density")
box on

