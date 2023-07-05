%% Script to randomly mix different populations of closed system crystallization and select 
% zircon from this probability distribution using the FC dZr and Bergell
clear all
clear rng
cd /Users/ltavazzani/Desktop/zurich_ID-TIMS/Keller_papers_zircSat/BayeZirChron.c-master/examples/BootstrappedCrystallizationDistributionExamples/UPB_LIT
% Load Bergell data and plot on rank-order plot
bergell=csvread('BergellAll.csv'); bergell(bergell==0)=NaN;
Fig1=figure(1); plot(bergell(:),'.','MarkerSize',10)
xlabel("N"); ylabel("Age, Ma")

% Load dZr data from FC simulation
dZr=[0.21,0.23,0.25,0.26,0.26,0.25,0.25,0.24,0.24,0.24,0.25,0.28,0.33,0.39,0.47,0.57,0.69,0.82,0.98,1.16,1.37,1.60,1.86,2.16,2.52,2.92,3.63,4.63,4.96,4.20,3.62,3.78,3.93,4.08,4.18,4.41,4.45,4.56,4.65,5.40,4.45,4.86,4.91,4.95,4.99,5.01,5.03,4.98,0.00];
dZr=[0.000000	0.000003	0.000003	0.000004	0.000004	0.000004	0.000004	0.000005	0.000004	0.000006	0.000006	0.000004	0.000004	0.000004	0.000002	0.000003	0.000002	0.000002	0.000001	0.000001	0.000002	0.000001	0.000001	0.000001	0.000000	0.000001	0.000001	0.000000	0.000001	0.000001	0.000000	0.000000	0.000001	0.000000	0.000001	0.000000	0.000000	0.000001	0.000000	0.000000	0.000000	0.000001	0.000000	0.000000	0.000000	0.000001	0.000000	0.000000	0.000000	0.000000	0.000001	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000001	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000001	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000];dZr=fliplr(dZr);
oL = length(dZr);

% Isolate Bergell samples and resample dZr for different time intervals 
bergell1=bergell(:,1); bergell1(isnan(bergell1)) = []; n=size(bergell1,1);
dZr1=interp1(1:oL, dZr, linspace(1,oL,n))';
bergell2=bergell(:,2); bergell2(isnan(bergell2)) = []; n=size(bergell2,1);
dZr2=interp1(1:oL, dZr, linspace(1,oL,n))';
bergell3=bergell(:,3);  bergell3(isnan(bergell3)) = []; n=size(bergell3,1);
dZr3=interp1(1:oL, dZr, linspace(1,oL,n))';
bergell4=bergell(:,4);  bergell4(isnan(bergell4)) = []; n=size(bergell4,1);
dZr4=interp1(1:oL, dZr, linspace(1,oL,n))';
bergell5=bergell(:,5); bergell5(isnan(bergell5)) = []; n=size(bergell5,1);
dZr5=interp1(1:oL, dZr, linspace(1,oL,n))';
bergell6=bergell(:,6);  bergell6(isnan(bergell6)) = []; n=size(bergell6,1);
dZr6=interp1(1:oL, dZr, linspace(1,oL,n))';
bergell7=bergell(:,7);  bergell7(isnan(bergell7)) = []; n=size(bergell7,1);
dZr7=interp1(1:oL, dZr, linspace(1,oL,n))';

% Plot resampled dZr for Bergell time intervals
Fig2=figure(2);
hold on;
plot(bergell1,dZr1); plot(bergell2,dZr2); plot(bergell3,dZr3);plot(bergell4,dZr4);
plot(bergell5,dZr5); plot(bergell6,dZr6); plot(bergell7,dZr7);
set(gca,'xdir','reverse')

% Addition of dZr in unique vector and random sampling of this
% vector using the function StoZirSamp

% Create cumulative matrix from all the resampled dZr
dZrC = nan(28,7);
dZrC(1:length(dZr1), 1) = dZr1; dZrC(1:length(dZr2), 2) = dZr2;
dZrC(1:length(dZr3), 3) = dZr3; dZrC(1:length(dZr4), 4) = dZr4;
dZrC(1:length(dZr5), 5) = dZr5; dZrC(1:length(dZr6), 6) = dZr6;
dZrC(1:length(dZr7), 7) = dZr7;

% Select proabilistic set-up
clear a A b B i ii kk x n N out P p t T TP TPORD Samp
Samp=3; % Number of random FC zircon crystallization intervals to sample
A=1; B=7;randomArray=(A-1)+(B-(A-1))*rand(1,Samp);
N=floor(randomArray)+1; N=unique(N,'first'); N=sort(N,'ascend');

% Extract random FC zircon crystallization intervals based on random sampling 
for ii=1:length(N)
t(:,ii) = bergell(:,N(ii));
end

for i=1:length(N)
p(:,i) = dZrC(:,N(i));
end

% Creating final randomly sampled matrices of time of crystallization (T)-probability
% of crystallization (P)
T=[]; P=[];
n=size(N,2);
for kk=1:length(n)
    if n==1
          T=vertcat(t(:,1));
          P=vertcat(p(:,1));
    elseif n==2
          T=vertcat(t(:,1),t(:,n(kk)));
          P=vertcat(p(:,1),p(:,n(kk)));
    elseif n==3
          T=vertcat(t(:,1),t(:,n(kk)),t(:,n(kk)-1));
          P=vertcat(p(:,1),p(:,n(kk)),p(:,n(kk)-1));
    elseif n==4
          T=vertcat(t(:,1),t(:,n(kk)),t(:,n(kk)-1),t(:,n(kk)-2)); 
          P=vertcat(p(:,1),p(:,n(kk)),p(:,n(kk)-1),p(:,n(kk)-2));
    elseif n==5
          T=vertcat(t(:,1),t(:,n(kk)),t(:,n(kk)-1),t(:,n(kk)-2),t(:,n(kk)-3));
          P=vertcat(p(:,1),p(:,n(kk)),p(:,n(kk)-1),p(:,n(kk)-2),p(:,n(kk)-3));
    elseif n==6
          T=vertcat(t(:,1),t(:,n(kk)),t(:,n(kk)-1),t(:,n(kk)-2),t(:,n(kk)-3),t(:,n(kk)-4));
          P=vertcat(p(:,1),p(:,n(kk)),p(:,n(kk)-1),p(:,n(kk)-2),p(:,n(kk)-3),p(:,n(kk)-4));
    else  
          T=vertcat(t(:,1),t(:,n(kk)),t(:,n(kk)-1),t(:,n(kk)-2),t(:,n(kk)-3),t(:,n(kk)-4),t(:,n(kk)-5));
          P=vertcat(p(:,1),p(:,n(kk)),p(:,n(kk)-1),p(:,n(kk)-2),p(:,n(kk)-3),p(:,n(kk)-4),p(:,n(kk)-5));
    end
T(isnan(T))=[];
P(isnan(P))=[];
end

TP=horzcat(T,P);TPORD=sortrows(TP,1,'descend');
Fig3=figure(3);
scatter(TPORD(:,1),TPORD(:,2),'r.');set(gca,'xdir','reverse')

%
clear t p n

% Time vector
t=TPORD(:,1);  
% Normalization of zircon crsytallization probabilities to 1
p=TPORD(:,2); p=smooth(p,20); 
p=p./sum(p);

%%%%%%%%%%%% Randomly generates a population of synthetic zircon crystallization ages
%%%%%%%%%%%% using the Matlab Built-in function randsrc()
nIt=100; nZr=30; alpha=0.1;
out = randsrc(nIt,nZr,[t';p']);

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
% Kernel Density Estimate of averaged zircon age spectra, scaled between 
% initiation and termination of zircon crystallization
[out_f,out_x,out_bw] = ksdensity(mean(out, 2));
outdistAvg = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,100));
x = linspace(0,1,100); outdistAvg = outdistAvg./trapz(x,outdistAvg);

% Plot average Kernel Density Estimate
Fig4=figure(4);
plot(x,outdistAvg,'Color','k','LineStyle','-','Linewidth',2)
hold on; 
% Kernel Density Estimates of all zircon age spectra, scaled between 
% initiation and termination of zircon crystallization
%for ii=1:nIt
%   [out_f,out_x,out_bw] = ksdensity(out(:,ii));
%   outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,100));
%   x = linspace(0,1,100); outdist = outdist./trapz(x,outdist); 
%   % Plot all Kernel Density Estimates
%   plot1=plot(x,outdist,'Color','k','LineStyle','-','Linewidth',0.15);
%   plot1.Color(4) = alpha;
%end

%% Compare to existing ID-TIMS datasets
cd /Users/ltavazzani/Desktop/LiteratureDatasets   
IDTIMS_DataLoad_KDE
% (**Volcanic eruptions**: Sesia_Caldera / Sesia_Caldera_NoAntecrysts  
%      Kneeling_Nun_Tuff / Huckleberry_Ridge_Tuff / Lava_Creek_Tuff 
%      Fish_Canyon_Tuff / Turkey_Creek_Caldera / Latir_Volcanics)
EmpDist_1=Kneeling_Nun_Tuff; 
EmpDist_2=Huckleberry_Ridge_Tuff;
EmpDist_3=Sesia_Caldera;
%EmpDist_4=Lava_Creek_Tuff;
% Insert igneous unit name for legend visualization (spelling does not matter)                       
EmpDist_1_str='Kneeling Nun Tuff'; 
EmpDist_2_str='Huckleberry Ridge Tuff'; 
EmpDist_3_str='Sesia Caldera';
%EmpDist_4_str='Lava Creek Tuff';
hold on; 
c = cool(4);
plot(x,EmpDist_1,'Color',c(1,:),'LineStyle','--','Linewidth',3)
plot(x,EmpDist_2,'Color',c(2,:),'LineStyle','-','Linewidth',3)
plot(x,EmpDist_3,'Color',c(3,:),'LineStyle','-.','Linewidth',3)
%plot(x,EmpDist_4,'Color',c(4,:),'LineStyle',':','Linewidth',3)
legend('Synthetic Zrn (avg.)',EmpDist_1_str,EmpDist_2_str,EmpDist_3_str,'AutoUpdate','off');
%legend('Synthetic Zrn (avg.)',EmpDist_1_str,EmpDist_2_str,EmpDist_3_str,EmpDist_4_str,'AutoUpdate','off');

% Plot details
ylim([0 3]);
set(gca,'xdir','reverse')
xlabel("Age (scaled)")
ylabel("Density")
box on
%close(Fig1,Fig2,Fig3)
%%
cumTPORD=(TPORD(:,2)./sum(TPORD(:,2))).*100
CUMP=cumsum(cumTPORD);
plot(TPORD(:,1),CUMP)
set(gca,'xdir','reverse')
%%
figure(8)
hold on;S=100;
plot(bergell1,smooth(dZr1,S)); 
plot(bergell2,smooth(dZr2,S)); 
plot(bergell3,smooth(dZr3,S));
plot(bergell4,smooth(dZr4,S)); 
plot(bergell5,smooth(dZr5,S));
plot(bergell6,smooth(dZr6,S)); 
plot(bergell7,smooth(dZr7,S));
set(gca,'xdir','reverse')
%%
tS9=[1.00	0.99	0.97	0.96	0.95	0.93	0.92	0.91	0.89	0.88	0.86	0.85	0.84	0.82	0.81	0.80	0.78	0.77	0.76	0.74	0.73	0.72	0.70	0.69	0.68	0.66	0.65	0.64	0.62	0.61	0.59	0.58	0.57	0.55	0.54	0.53	0.51	0.50	0.49	0.47	0.46	0.45	0.43	0.42	0.41	0.39	0.38	0.36	0.35	0.34	0.32	0.31	0.30	0.28	0.27	0.26	0.24	0.23	0.22	0.20	0.19	0.18	0.16	0.15	0.14	0.12	0.11	0.09	0.08	0.07	0.05	0.04	0.03	0.01	0.00];
dZrS9=[0.0037180	0.0115560	0.0111971	0.0108490	0.0105114	0.0101840	0.0098665	0.0095585	0.0092597	0.0084576	0.0103427	0.0084196	0.0081556	0.0078994	0.0081376	0.0073817	0.0072245	0.0069939	0.0067702	0.0065532	0.0000000	0.0072761	0.0071665	0.0000000	0.0000000	0.0085877	0.0083918	0.0081999	0.0080119	0.0078278	0.0076475	0.0272723	0.0133410	0.0166769	0.0107378	0.0095293	0.0065671	0.0110957	0.0088172	0.0085357	0.0082625	0.0000000	0.0081369	0.0080131	0.0000000	0.0072005	0.0095401	0.0093205	0.0091054	0.0088948	0.0141617	0.0339322	0.0172608	0.0148157	0.0111383	0.0077287	0.0128390	0.0103019	0.0099726	0.0000000	0.0090096	0.0088714	0.0000000	0.0098566	0.0105079	0.0102644	0.0100259	0.0204237	0.0358252	0.0232003	0.0143939	0.0097306	0.0140391	0.0119827	0.0115998];
hold on;S=0.06;
plot(tS9,smooth(dZrS9,S)); 
set(gca,'xdir','reverse')
