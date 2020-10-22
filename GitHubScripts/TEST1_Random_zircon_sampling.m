%% Randomly generated zircon samples based on MCS outputs of relative zircon 
% growth under different crystallization scenarios 
% 1/7/2020 Lorenzo Tavazzani

clear all
close all
clc

% Data from literature used in bootstrapped crystallization dist. plots
cd /Users/Utente/Desktop/zurich_ID-TIMS/Keller_papers_zircSat/BayeZirChron.c-master/examples/BootstrappedCrystallizationDistributionExamples/SCRIPTS
IDTIMS_DataLoad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sample for statistical comparison

nSimulations=100; 
nZircons=10000;

%(TEST2fc - Mode 3)
t=[22 21 20	19 18 17 16	15 14 13 12	11 10 9 8 7 6 5 4]; % population
p=[9.99	12.48 12.27	12.18 11.33	10.61 9.74 10.74 9.71 7.06 4.49	3.03 1.97 1.23 0.79	0.61 0.61 0.64 0.61]; % probabilities
p=smooth(p)'; % smoothing filter on the data, can be avoided
p=p./sum(p); % normalization of growth rates to p=1

% randomly generated population
out = randsrc(nSimulations,nZircons,[t;p]);

% order "youngest" to "oldest"
a=size(out); a=a(:,1);
for ii=1:a
out(ii,:)=sort(out(ii,:),'ascend');
end
out=out';

for i=1:size(out,2)
    out(:,i) = out(:,i) - min(out(:,i));
    out(:,i) = out(:,i)./max(out(:,i));
end

% Bootstrapped 
npoints = 100;
[out_f,out_x,out_bw] = ksdensity(mean(out, 2));
outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,npoints));
x = linspace(0,1,npoints);
outdistKSTest = outdist./trapz(x,outdist);

%(TEST2a - Mode 2)
tFC1=[7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	33	34	35	36	37	38	39]; % population
pFC1=[1	1	1	1	1	3	4	0	3	6	14	41	17	8	11	12	13	15	19	11	7	8	9	9	10	11	4]; % probabilities
pFC1=smooth(pFC1)'; % smoothing filter on the data, can be avoided
pFC1=pFC1./sum(pFC1); % normalization of growth rates to p=1

% randomly generated population
outFC1 = randsrc(nSimulations,nZircons,[tFC1;pFC1]);

% order "youngest" to "oldest"
a=size(outFC1); a=a(:,1);
for ii=1:a
outFC1(ii,:)=sort(outFC1(ii,:),'ascend');
end
outFC1=outFC1';

for i=1:size(outFC1,2)
    outFC1(:,i) = outFC1(:,i) - min(outFC1(:,i));
    outFC1(:,i) = outFC1(:,i)./max(outFC1(:,i));
end

% Bootstrapped 
npoints = 100;
[outFC1_f,outFC1_x,outFC1_bw] = ksdensity(mean(outFC1, 2));
outFC1dist = interp1(outFC1_x,outFC1_f,linspace(0-0.05,1+outFC1_bw,npoints));
x = linspace(0,1,npoints);
outdistKSTestFC1 = outFC1dist./trapz(x,outFC1dist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Randomly generated samples for different modeling conditions
close all 
% Randomly generated samples for FC results 
% (TEST2fc - Mode 3)
% on sample PST214 composition

nSimulations=100; 
alpha=0.1; % Adjust transparency based on number of simulations (nSimulations)
nZircons=94;

%Complete crystallization
t=[22 21 20	19 18 17 16	15 14 13 12	11 10 9 8 7 6 5 4]; % population
p=[9.99	12.48 12.27	12.18 11.33	10.61 9.74 10.74 9.71 7.06 4.49	3.03 1.97 1.23 0.79	0.61 0.61 0.64 0.61]; % probabilities
p=smooth(p)'; % smoothing filter on the data, can be avoided
p=p./sum(p); % normalization of growth rates to p=1

% randomly generated population
out = randsrc(nSimulations,nZircons,[t;p]);

% order "youngest" to "oldest"
a=size(out); a=a(:,1);
for ii=1:a
out(ii,:)=sort(out(ii,:),'ascend');
end
out=out';
figure(1)

% histogram of randomly generated samples
subplot(2,2,1)
b=size(t,2); 
hold on
for i=1:size(out,2)
   histogram(out(:,i),linspace(min(t),max(t),b)); hold on
end

xlim([min(t)-1 max(t)+1]); ylim([0 (nZircons/2)])
set(gca,'xdir','reverse'); box on
xlabel("t interval"); ylabel("N"); 

% Plot (rank order) data
subplot(2,2,3)
plot(out(:),'.','MarkerSize',10)
ylim([min(t) max(t)])
xlabel("N"); ylabel("Age, Ma")

% Plot (rank order) scaled data
subplot(2,2,4)
for i=1:size(out,2)
    out(:,i) = out(:,i) - min(out(:,i));
    out(:,i) = out(:,i)./max(out(:,i));
end
plot(out(:),'.','MarkerSize',10)
xlabel("N"); ylabel("Age, scaled")

% Bootstrapped 
subplot(2,2,2)
hold on    
c = lines(7);
plot(x,bdist,'Color',c(1,:),'LineStyle','-','Linewidth',2)
%plot(x,edist,'Color',c(6,:),'LineStyle','-.','Linewidth',2)
%plot(x,Sdist,'Color',c(3,:),'LineStyle','-.','Linewidth',1)

%MELTS_Keller,2018
%Nt=[1.00	0.99	0.99	0.98	0.98	0.97	0.96	0.95	0.95	0.94	0.94	0.93	0.92	0.91	0.90	0.89	0.88	0.87	0.85	0.84	0.82	0.80	0.79	0.77	0.76	0.75	0.73	0.72	0.70	0.69	0.68	0.67	0.65	0.64	0.63	0.62	0.60	0.59	0.57	0.56	0.54	0.53	0.52	0.50	0.48	0.47	0.45	0.43	0.41	0.39	0.37	0.36	0.35	0.33	0.31	0.29	0.28	0.26	0.24	0.22	0.19	0.17	0.15	0.13	0.12	0.09	0.07	0.05	0.03	0.00];
%NMZr=[0.03	0.18	0.34	0.45	0.57	0.74	0.92	1.08	1.25	1.39	1.53	1.70	1.85	2.05	2.19	2.32	2.39	2.40	2.40	2.34	2.26	2.18	2.09	2.02	1.94	1.87	1.79	1.71	1.62	1.55	1.49	1.45	1.39	1.34	1.29	1.24	1.19	1.14	1.08	1.03	0.99	0.94	0.89	0.86	0.82	0.78	0.73	0.69	0.64	0.61	0.58	0.56	0.53	0.51	0.49	0.46	0.44	0.42	0.40	0.37	0.34	0.32	0.30	0.28	0.27	0.25	0.23	0.22	0.21	0.19];
%plot(smooth(Nt),smooth(NMZr),'Color','b','LineStyle','-','Linewidth',1)

%Bergell_Keller,2018
%Nt=[1.00	0.99	0.99	0.98	0.98	0.97	0.97	0.96	0.96	0.95	0.95	0.94	0.94	0.93	0.93	0.92	0.92	0.91	0.90	0.90	0.89	0.89	0.88	0.87	0.87	0.86	0.85	0.85	0.84	0.83	0.82	0.81	0.80	0.79	0.78	0.77	0.76	0.74	0.73	0.72	0.71	0.70	0.69	0.68	0.67	0.66	0.65	0.64	0.63	0.63	0.62	0.61	0.60	0.59	0.58	0.57	0.56	0.55	0.54	0.53	0.51	0.50	0.49	0.48	0.47	0.45	0.44	0.43	0.42	0.40	0.39	0.37	0.36	0.34	0.33	0.32	0.30	0.28	0.26	0.24	0.22	0.20	0.19	0.17	0.15	0.13	0.11	0.08	0.06	0.04	0.02	0.00];
%NMZr=[0.49	0.51	0.54	0.57	0.61	0.66	0.70	0.74	0.78	0.83	0.88	0.92	0.98	1.01	1.05	1.10	1.15	1.20	1.25	1.30	1.33	1.38	1.43	1.47	1.52	1.56	1.60	1.65	1.69	1.72	1.76	1.79	1.81	1.83	1.85	1.85	1.85	1.85	1.84	1.82	1.80	1.78	1.75	1.72	1.69	1.67	1.64	1.61	1.58	1.55	1.51	1.48	1.46	1.43	1.40	1.36	1.32	1.29	1.24	1.20	1.14	1.10	1.06	1.03	0.99	0.96	0.92	0.90	0.85	0.81	0.78	0.74	0.71	0.68	0.64	0.61	0.58	0.54	0.50	0.46	0.44	0.42	0.41	0.40	0.39	0.38	0.36	0.33	0.30	0.26	0.22	0.17];
%plot(smooth(Nt),smooth(NMZr),'Color','m','LineStyle','-','Linewidth',1)

%Watson_Keller,2018
Nt=[1.00	1.00	0.99	0.99	0.99	0.98	0.98	0.98	0.97	0.97	0.96	0.96	0.95	0.95	0.93	0.92	0.91	0.90	0.88	0.87	0.86	0.84	0.82	0.81	0.79	0.78	0.76	0.75	0.74	0.73	0.71	0.70	0.69	0.67	0.66	0.65	0.63	0.62	0.60	0.58	0.56	0.55	0.54	0.52	0.51	0.49	0.47	0.46	0.44	0.43	0.41	0.40	0.39	0.37	0.35	0.34	0.33	0.31	0.30	0.28	0.26	0.24	0.22	0.20	0.19	0.17	0.14	0.13	0.11	0.09	0.08	0.06	0.04	0.03	0.02	0.00];
NMZr=[0.32	0.41	0.53	0.64	0.74	0.83	0.92	0.99	1.07	1.16	1.27	1.35	1.45	1.50	1.58	1.65	1.72	1.76	1.80	1.83	1.86	1.86	1.85	1.83	1.81	1.78	1.75	1.73	1.70	1.66	1.64	1.60	1.57	1.53	1.49	1.45	1.41	1.36	1.31	1.26	1.21	1.16	1.12	1.06	1.03	0.99	0.94	0.90	0.87	0.83	0.79	0.76	0.73	0.69	0.66	0.63	0.61	0.58	0.56	0.52	0.50	0.46	0.43	0.41	0.38	0.36	0.34	0.32	0.29	0.28	0.27	0.25	0.24	0.22	0.21	0.20];
plot(smooth(Nt),smooth(NMZr),'Color','r','LineStyle','-','Linewidth',2)

npoints = 100;
[out_f,out_x,out_bw] = ksdensity(mean(out, 2));
outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,npoints));
x = linspace(0,1,npoints);
outdist = outdist./trapz(x,outdist);
plot(x,outdist,'Color','k','LineStyle','-','Linewidth',2)

%[h,p] = kstest2(outdist,bdist,'Alpha',0.02) %Kolmolarov-Smirnov test

hold on; 
for ii=1:nSimulations
   [out_f,out_x,out_bw] = ksdensity(out(:,ii));
   outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,npoints));
   x = linspace(0,1,npoints);
   outdist = outdist./trapz(x,outdist); 
   [h,p] = kstest2(outdist,outdistKSTest,'Alpha',0.01); %Kolmolarov-Smirnov test
   Tresult(ii) = h;
   Presult(ii) = p; Presult=sort(Presult,'descend');
   plot1=plot(x,outdist,'Color','k','LineStyle','-','Linewidth',0.15);
   plot1.Color(4) = alpha;
end

PTestFC=size(find(Tresult==0),2)

legend('Samperton - Bergell','Watson (1996)','Synthetic Zrn');

ylim([0 3]);
set(gca,'xdir','reverse')
xlabel("Age (scaled)")
ylabel("Density")
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Randomly generated samples for FC with 1 silicic recharges results
% (TEST2a - Mode 2)
% on sample PST214 composition

%Complete crystallization
t=[7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	33	34	35	36	37	38	39]; % population
p=[1	1	1	1	1	3	4	0	3	6	14	41	17	8	11	12	13	15	19	11	7	8	9	9	10	11	4]; % probabilities
%Cutoff at 750 °C (eruption)
t=[20	21	22	23	24	25	26	33	34	35	36	37	38	39]; % population
p=[8	11	12	13	15	19	11	7	8	9	9	10	11	4]; % probabilities
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
subplot(2,2,1)
b=size(t,2); 
hold on
for i=1:size(out,2)
   histogram(out(:,i),linspace(min(t),max(t),max(t)-min(t)+1)); hold on
end

xlim([min(t)-1 max(t)+1]); ylim([0 nZircons/2])
set(gca,'xdir','reverse'); box on
xlabel("t interval"); ylabel("N"); 

% % Plot (rank order) data
subplot(2,2,3)
plot(out(:),'.','MarkerSize',10)
ylim([min(t) max(t)])
xlabel("N"); ylabel("Age (No unit)")

% Plot (rank order) scaled data
subplot(2,2,4)
for i=1:size(out,2)
    out(:,i) = out(:,i) - min(out(:,i));
    out(:,i) = out(:,i)./max(out(:,i));
end
plot(out(:),'.','MarkerSize',10)
xlabel("N"); ylabel("Age, scaled")

%Bootstrapped 
subplot(2,2,2)

hold on; 
plot(x,Sdist,'Color',c(2,:),'LineStyle','--','Linewidth',1)
plot(x,Kdist,'Color',c(3,:),'LineStyle','-.','Linewidth',1)
plot(x,HRTdist,'Color',c(7,:),'LineStyle',':','Linewidth',1)

npoints = 100;
[out_f,out_x,out_bw] = ksdensity(mean(out, 2));
outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,npoints));
x = linspace(0,1,npoints);
outdist = outdist./trapz(x,outdist);
plot(x,outdist,'Color','k','LineStyle','-','Linewidth',2)

hold on; 
for ii=1:nSimulations
   [out_f,out_x,out_bw] = ksdensity(out(:,ii));
   outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,npoints));
   x = linspace(0,1,npoints);
   outdist = outdist./trapz(x,outdist); 
   [h,p] = kstest2(outdist,outdistKSTestFC1,'Alpha',0.01); %Kolmolarov-Smirnov test
   Tresult(ii) = h;
   Presult(ii) = p; Presult=sort(Presult,'descend');
   plot1=plot(x,outdist,'Color','k','LineStyle','-','Linewidth',0.15);
   plot1.Color(4) = alpha;
end

PTestFC1=size(find(Tresult==0),2)

legend('SesiaG','KNT','HRT','Synthetic Zrn','Location','northwest');

ylim([0 3]);
set(gca,'xdir','reverse')
xlabel("Age (scaled)")
ylabel("Density")
box on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Randomly generated samples for FC with 2 silicic recharges results
% (TEST2c - Mode 2)
% on sample PST214 composition

%Complete crystallization
t=[7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	31	32	33	34	35	36	43	44	45	46	47	48	49]; % population
p=[1	1	1	1	2	2	2	3	4	8	11	42	38	13	17	14	17	25	23	15	2	11	12	13	15	19	11	7	8	9	9	10	11	4]; % probabilities
%Cutoff at 750 °C (eruption)
t=[	21	22	23	24	25	26	27	31	32	33	34	35	36	43	44	45	46	47	48	49]; % population
p=[	17	14	17	25	23	15	2	11	12	13	15	19	11	7	8	9	9	10	11	4]; % probabilities

p=smooth(p)'; % smoothing filter on the data, can be avoided
p=p./sum(p);  % normalization of growth rates to p=1

% randomly generated population
out = randsrc(nSimulations,nZircons,[t;p]);

% order "youngest" to "oldest"
a=size(out); a=a(:,1);
for ii=1:a
out(ii,:)=sort(out(ii,:),'ascend');
end

figure(3)
% histogram of randomly generated samples
out=out';
subplot(2,2,1)
b=size(t,2); 
hold on
for i=1:size(out,2)
   histogram(out(:,i),linspace(min(t),max(t),max(t)-min(t)+1)); hold on
end

xlim([min(t)-1 max(t)+1]); ylim([0 nZircons/2])
set(gca,'xdir','reverse'); box on
xlabel("t interval"); ylabel("N"); 

% % Plot (rank order) data
subplot(2,2,3)
plot(out(:),'.','MarkerSize',10)
ylim([min(t) max(t)])
xlabel("N"); ylabel("Age (No unit)")

% Plot (rank order) scaled data
subplot(2,2,4)
for i=1:size(out,2)
    out(:,i) = out(:,i) - min(out(:,i));
    out(:,i) = out(:,i)./max(out(:,i));
end
plot(out(:),'.','MarkerSize',10)
xlabel("N"); ylabel("Age, scaled")

%Bootstrapped 
subplot(2,2,2)

hold on; 
plot(x,Sdist,'Color',c(2,:),'LineStyle','--','Linewidth',1)
plot(x,Kdist,'Color',c(3,:),'LineStyle','-.','Linewidth',1)
plot(x,HRTdist,'Color',c(7,:),'LineStyle',':','Linewidth',1)

npoints = 100;
[out_f,out_x,out_bw] = ksdensity(mean(out, 2));
outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,npoints));
x = linspace(0,1,npoints);
outdist = outdist./trapz(x,outdist);
plot(x,outdist,'Color','k','LineStyle','-','Linewidth',2)

hold on; 
for ii=1:nSimulations
   [out_f,out_x,out_bw] = ksdensity(out(:,ii));
   outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,npoints));
   x = linspace(0,1,npoints);
   outdist = outdist./trapz(x,outdist); 
   plot1=plot(x,outdist,'Color','k','LineStyle','-','Linewidth',0.15);
   plot1.Color(4) = alpha;
end

legend('SesiaG','KNT','HRT','Synthetic Zrn','Location','northwest');

ylim([0 3]);
set(gca,'xdir','reverse')
xlabel("Age (scaled)")
ylabel("Density")
box on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Randomly generated samples for FC with 3 silicic recharges results
% (TEST2d - Mode 2)
% on sample PST214 composition

%Complete crystallization
t=[8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	30	31	32	33	34	35	36	40	41	42	43	44	45	52	53	54	55	56	57	58]; % population
p=[1	2	2	2	2	3	4	5	7	14	37	62	17	18	18	29	32	17	17	0	16	17	14	17	25	23	15	2	11	12	13	15	19	11	7	8	9	9	10	11	4]; % probabilities
%Cutoff at 750 °C (eruption)
t=[	21	22	23	24	25	26	27	28	30	31	32	33	34	35	36	40	41	42	43	44	45	52	53	54	55	56	57	58]; % population
p=[	18	18	29	32	17	17	0	16	17	14	17	25	23	15	2	11	12	13	15	19	11	7	8	9	9	10	11	4]; % probabilities
p=smooth(p)'; % smoothing filter on the data, can be avoided
p=p./sum(p);  % normalization of growth rates to p=1

% randomly generated population
out = randsrc(nSimulations,nZircons,[t;p]);

% order "youngest" to "oldest"
a=size(out); a=a(:,1);
for ii=1:a
out(ii,:)=sort(out(ii,:),'ascend');
end

figure(4)
% histogram of randomly generated samples
out=out';
subplot(2,2,1)
b=size(t,2); 
hold on
for i=1:size(out,2)
   histogram(out(:,i),linspace(min(t),max(t),max(t)-min(t)+1)); hold on
end

xlim([min(t)-1 max(t)+1]); ylim([0 nZircons/2])
set(gca,'xdir','reverse'); box on
xlabel("t interval"); ylabel("N"); 

% Plot (rank order) data
subplot(2,2,3)
plot(out(:),'.','MarkerSize',10)
ylim([min(t) max(t)])
xlabel("N"); ylabel("Age (No unit)")

% Plot (rank order) scaled data
subplot(2,2,4)
for i=1:size(out,2)
    out(:,i) = out(:,i) - min(out(:,i));
    out(:,i) = out(:,i)./max(out(:,i));
end
plot(out(:),'.','MarkerSize',10)
xlabel("N"); ylabel("Age, scaled")


%Bootstrapped 
subplot(2,2,2)

hold on; 
plot(x,Sdist,'Color',c(2,:),'LineStyle','--','Linewidth',1)
%plot(x,Kdist,'Color',c(3,:),'LineStyle','-.','Linewidth',1)
%plot(x,HRTdist,'Color',c(7,:),'LineStyle',':','Linewidth',1)

npoints = 100;
[out_f,out_x,out_bw] = ksdensity(mean(out, 2));
outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,npoints));
x = linspace(0,1,npoints);
outdist = outdist./trapz(x,outdist);
plot(x,outdist,'Color','k','LineStyle','-','Linewidth',2)

hold on; 
for ii=1:nSimulations
   [out_f,out_x,out_bw] = ksdensity(out(:,ii));
   outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,npoints));
   x = linspace(0,1,npoints);
   outdist = outdist./trapz(x,outdist); 
   plot1=plot(x,outdist,'Color','k','LineStyle','-','Linewidth',0.15);
   plot1.Color(4) = alpha;
end

legend('SesiaG','KNT','HRT','Synthetic Zrn','Location','northwest');

ylim([0 3]);
set(gca,'xdir','reverse')
xlabel("Age (scaled)")
ylabel("Density")
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Randomly generated samples for FC with 2 silicic and 1 mafic recharges results
% (TEST2e - Mode 2)
% on sample PST214 composition

t=[8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	29	30	31	32	33	34	35	39	40	41	42	43	44	51	52	53	54	55	56	57]; % population
p=[2	2	2	3	7	7	39	104	13	14	12	12	13	14	10	17	14	17	25	23	15	2	11	12	13	15	19	11	7	8	9	9	10	11	4]; % probabilities
%Cutoff at 750 °C (eruption)
%t=[	21	22	29	30	31	32	33	34	35	39	40	41	42	43	44	51	52	53	54	55	56	57]; % population
%p=[	14	10	17	14	17	25	23	15	2	11	12	13	15	19	11	7	8	9	9	10	11	4]; % probabilities
p=smooth(p)'; % smoothing filter on the data, can be avoided
p=p./sum(p);  % normalization of growth rates to p=1

% randomly generated population
out = randsrc(nSimulations,nZircons,[t;p]);

% order "youngest" to "oldest"
a=size(out); a=a(:,1);
for ii=1:a
out(ii,:)=sort(out(ii,:),'ascend');
end

figure(5)
% histogram of randomly generated samples
out=out';
subplot(2,2,1)
b=size(t,2); 
hold on
for i=1:size(out,2)
   histogram(out(:,i),linspace(min(t),max(t),max(t)-min(t)+1)); hold on
end

xlim([min(t)-1 max(t)+1]); ylim([0 nZircons/2])
set(gca,'xdir','reverse'); box on
xlabel("t interval"); ylabel("N"); 

% Plot (rank order) data
subplot(2,2,3)
plot(out(:),'.','MarkerSize',10)
ylim([min(t) max(t)])
xlabel("N"); ylabel("Age (No unit)")

% Plot (rank order) scaled data
subplot(2,2,4)
for i=1:size(out,2)
    out(:,i) = out(:,i) - min(out(:,i));
    out(:,i) = out(:,i)./max(out(:,i));
end
plot(out(:),'.','MarkerSize',10)
xlabel("N"); ylabel("Age, scaled")

%Bootstrapped 
subplot(2,2,2)

hold on; 
plot(x,Sdist,'Color',c(2,:),'LineStyle','--','Linewidth',1)
plot(x,Kdist,'Color',c(3,:),'LineStyle','-.','Linewidth',1)
plot(x,HRTdist,'Color',c(7,:),'LineStyle',':','Linewidth',1)

npoints = 100;
[out_f,out_x,out_bw] = ksdensity(mean(out, 2));
outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,npoints));
x = linspace(0,1,npoints);
outdist = outdist./trapz(x,outdist);
plot(x,outdist,'Color','k','LineStyle','-','Linewidth',2)

hold on; 
for ii=1:nSimulations
   [out_f,out_x,out_bw] = ksdensity(out(:,ii));
   outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,npoints));
   x = linspace(0,1,npoints);
   outdist = outdist./trapz(x,outdist); 
   plot1=plot(x,outdist,'Color','k','LineStyle','-','Linewidth',0.15);
   plot1.Color(4) = alpha;
end

legend('SesiaG','KNT','HRT','Synthetic Zrn','Location','northwest');

ylim([0 3]);
set(gca,'xdir','reverse')
xlabel("Age (scaled)")
ylabel("Density")
box on