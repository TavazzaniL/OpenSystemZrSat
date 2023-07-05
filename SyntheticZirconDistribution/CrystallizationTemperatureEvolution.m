%%% Random sampling of zircon crystallization temperatures 

cd /Users/ltavazzani/Desktop/SyntheticZirconDistribution

                         %%%%%% [User Input] %%%%%%  
sim=14; mod=0; nIt=1; nZr=30; cycles=3; alpha = 0.2; G=25;%+-uncertainity(°C)

%%% Random selection of zircon crystallization steps and temperatures
Out=ones(nZr,cycles); 
Out1=ones(nZr,cycles);
for c = 1:cycles 
    [out,out1] = StoZrnSamp_Temp(sim,mod,nIt,nZr);
    Out(:,c)=out;
    Out1(:,c)=out1;
end

%%% Addition of random error
G=50*G^-1;
g=round(((rand(nZr,cycles)-1/2)*100)./G); 
Out1=Out1+g;

O=Out(:); O1=Out1(:);
OO=horzcat(O,O1); OO=sortrows(OO);


%%%%%% Temperature evolution plot

Cent=[0.95 0.85 0.75 0.65 0.55 0.45 0.35 0.25 0.15 0.05]; %Centroids diagrams
box1=find(OO(:,1)<=1.0 & OO(:,1)>0.9);
box2=find(OO(:,1)<0.9 & OO(:,1)>0.8);
box3=find(OO(:,1)<0.8 & OO(:,1)>0.7);
box4=find(OO(:,1)<0.7 & OO(:,1)>0.6); 
box5=find(OO(:,1)<0.6 & OO(:,1)>0.5); 
box6=find(OO(:,1)<0.5 & OO(:,1)>0.4); 
box7=find(OO(:,1)<0.4 & OO(:,1)>0.3); 
box8=find(OO(:,1)<0.3 & OO(:,1)>0.2); 
box9=find(OO(:,1)<0.2 & OO(:,1)>0.1); 
box10=find(OO(:,1)<0.1 & OO(:,1)>=0.0); 

dc = [0 0.4470 0.7410];
dr = [0.8500, 0.3250, 0.0980];
if sim==1 || sim==23
    cc=dc;
else
    cc=dr;
end

%%% Box and whisker plot to summarize temperature evolution
hFig = figure(1);
subplot(3,2,1)
set(hFig, 'Position', [200 200 824 486])
El=OO(:,2); ElLabel='T (°C)'; ElLim=[650 850]; 
C = [El(box1)' El(box2)' El(box3)' El(box4)' El(box5)' El(box6)' El(box7)' El(box8)' El(box9)' El(box10)'];
grp = [Cent(1).*ones(1,size(box1,1)),Cent(2).*ones(1,size(box2,1)),Cent(3).*ones(1,size(box3,1)),Cent(4).*ones(1,size(box4,1)),Cent(5).*ones(1,size(box5,1)),Cent(6).*ones(1,size(box6,1)),Cent(7).*ones(1,size(box7,1)),Cent(8).*ones(1,size(box8,1)),Cent(9).*ones(1,size(box9,1)),Cent(10).*ones(1,size(box10,1))];
boxplot(C,grp); 
set(gca,'xdir','reverse'); xlabel("t interval"); set(gca,'ylim',ElLim);

%%% Temperature values
subplot(3,2,3)
scatter(OO(:,1),OO(:,2),'MarkerEdgeColor',cc)   
ylim([650 850])
set(gca,'xdir','reverse'); box on
xlabel("t interval"); ylabel("T(°C)");

%%% Histogram and Kernel density of temperature values
subplot(3,2,4)
bw=5;
edges=[650:5:850];
[counts,bins]=hist(OO(:,2),edges); 
barh(bins,counts,'FaceColor',cc,'HandleVisibility','off'); hold on;
[f,xi] = ksdensity(OO(:,2),'bandwidth',bw);
plot(f*500,xi,'k-','Color',cc,'Linewidth',2);hold on;
legend('All')
ylim([650 850]); xlim([0 20])
box on; %set(gca,'xdir','reverse')
xlabel("N");ylabel("T(°C)")
                 
%%%%%% KDE plot
% Kernel Density Estimate of averaged zircon age spectra, scaled between 
% initiation and termination of zircon crystallization
out=Out;
nIt = cycles; 
[out_f,out_x,out_bw] = ksdensity(mean(out, 2));
outdistAvg = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,100));
x = linspace(0,1,100); outdistAvg = outdistAvg./trapz(x,outdistAvg);

% Plot average Kernel Density Estimate
subplot(3,2,5)
plot(x,outdistAvg,'LineStyle','-','Linewidth',1)
hold on; 
legend('Synthetic Zrn (avg.)','AutoUpdate','off');

% Kernel Density Estimates of all zircon age spectra, scaled between 
% initiation and termination of zircon crystallization
for ii=1:nIt
   [out_f,out_x,out_bw] = ksdensity(out(:,ii));
   outdist = interp1(out_x,out_f,linspace(0-0.05,1+out_bw,100));
   x = linspace(0,1,100); outdist = outdist./trapz(x,outdist); 
   % Plot all Kernel Density Estimates
   plot1=plot(x,outdist,'Color',cc,'LineStyle','-','Linewidth',0.15);
   plot1.Color(4) = alpha;
end

% Plot details
ylim([0 3]);
set(gca,'xdir','reverse')
xlabel("Age (scaled)")
ylabel("Density")
box on
