%% Load and calculate Kernel Density Estimates (KDEs) of empirical zircon age 
% distributions from Tavazzani et al., 2020 and literature 
% 
% Author: Lorenzo Tavazzani
% e-mail address: ltavazzani@smu.edu
% Release: 1.0
% Release date: 10/23/20                  


                          %%%%%% [User Input] %%%%%% 
% Select path to folder that contain the scripts and .csv files of 
% data from literature 
cd /Users/Utente/Desktop/LiteratureDatasets                       


%%%%%%% [1] Sesia Magmatic System (Karakas et al., 2019; Tavazzani et al., 2020)
% Valle Mosso pluton
SesiaG=csvread('SesiaGranites.csv');
SesiaG(SesiaG==0)=NaN;
for i=1:size(SesiaG,2)
    SesiaG(:,i) = SesiaG(:,i) - min(SesiaG(:,i));
    SesiaG(:,i) = SesiaG(:,i)./max(SesiaG(:,i));
end
% Valle Mosso pluton (antecrysts excluded)
SesiaGNA=csvread('SesiaGranitesNoANT.csv');
SesiaGNA(SesiaGNA==0)=NaN;
for i=1:size(SesiaGNA,2)
    SesiaGNA(:,i) = SesiaGNA(:,i) - min(SesiaGNA(:,i));
    SesiaGNA(:,i) = SesiaGNA(:,i)./max(SesiaGNA(:,i));
end
% Sesia Caldera rhyolites
SesiaV=csvread('SesiaVolcanics.csv');
SesiaV(SesiaV==0)=NaN;
for i=1:size(SesiaV,2)
    SesiaV(:,i) = SesiaV(:,i) - min(SesiaV(:,i));
    SesiaV(:,i) = SesiaV(:,i)./max(SesiaV(:,i));
end
% Sesia Caldera rhyolites (antecrysts excluded)
SesiaVNA=csvread('SesiaVolcanicsNoANT.csv');
SesiaVNA(SesiaVNA==0)=NaN;
for i=1:size(SesiaVNA,2)
    SesiaVNA(:,i) = SesiaVNA(:,i) - min(SesiaVNA(:,i));
    SesiaVNA(:,i) = SesiaVNA(:,i)./max(SesiaVNA(:,i));
end

%%%%%%% [2] Bergell pluton (Samperton et al., 2015)
bergell=csvread('BergellAll.csv');
bergell(bergell==0)=NaN;
for i=1:size(bergell,2)
    bergell(:,i) = bergell(:,i) - min(bergell(:,i));
    bergell(:,i) = bergell(:,i)./max(bergell(:,i));
end

%%%%%%% [3] Mt Capanne pluton (Barboni et al., 2015)
elba = csvread('ElbaAll.csv');
elba(elba==0)=NaN;
for i=1:size(elba,2)
    elba(:,i) = elba(:,i) - min(elba(:,i));
    elba(:,i) = elba(:,i)./max(elba(:,i));
end

%%%%%%% [4] Mt Givens granodiorite (Frazer et al., 2014)
MtGivens = csvread('MtGivensAll.csv');
MtGivens(MtGivens==0)=NaN;
for i=1:size(MtGivens,2)
    MtGivens(:,i) = MtGivens(:,i) - min(MtGivens(:,i));
    MtGivens(:,i) = MtGivens(:,i)./max(MtGivens(:,i));
end

%%%%%%% [5] Mt Stuart Intrusive complex and Tenpeak pluton (Matzel et al., 2007)
% Mt Stuart Intrusive complex
MtStuart = csvread('MtStuartAll.csv');
MtStuart(MtStuart==0)=NaN;
for i=1:size(MtStuart,2)
    MtStuart(:,i) = MtStuart(:,i) - min(MtStuart(:,i));
    MtStuart(:,i) = MtStuart(:,i)./max(MtStuart(:,i));
end
% Tenpeak pluton
TPeak = csvread('TenpeakAll.csv');
TPeak(TPeak==0)=NaN;
for i=1:size(TPeak,2)
    TPeak(:,i) = TPeak(:,i) - min(TPeak(:,i));
    TPeak(:,i) = TPeak(:,i)./max(TPeak(:,i));
end

%%%%%%% [6] Tuolomne Intrusive suite (Coleman et al., 2004)
Tuolomne = csvread('TuolomneAll.csv');
Tuolomne(Tuolomne==0)=NaN;
for i=1:size(Tuolomne,2)
    Tuolomne(:,i) = Tuolomne(:,i) - min(Tuolomne(:,i));
    Tuolomne(:,i) = Tuolomne(:,i)./max(Tuolomne(:,i));
end

%%%%%%% [7] John Muir Intrusive suite (Davis et al., 2012)
JMuir = csvread('JohnMuirAll.csv');
JMuir(JMuir==0)=NaN;
for i=1:size(JMuir,2)
    JMuir(:,i) = JMuir(:,i) - min(JMuir(:,i));
    JMuir(:,i) = JMuir(:,i)./max(JMuir(:,i));
end

%%%%%%% Southern Adamello Batholith 
% [8] Lago della Vacca (Schoene et al., 2012)
LDVC = csvread('LagoVaccaAll.csv');
LDVC(LDVC==0)=NaN;
for i=1:size(LDVC,2)
    LDVC(:,i) = LDVC(:,i) - min(LDVC(:,i));
    LDVC(:,i) = LDVC(:,i)./max(LDVC(:,i));
end
% [9] Val Fredda Igneous Complex (Broderick et al., 2015)
VFIC = csvread('ValFreddaAll.csv');
VFIC(VFIC==0)=NaN;
for i=1:size(VFIC,2)
    VFIC(:,i) = VFIC(:,i) - min(VFIC(:,i));
    VFIC(:,i) = VFIC(:,i)./max(VFIC(:,i));
end

%%%%%%% [10] La Gloria pluton (Gutierrez et al. 2018)
lagloria=csvread('LaGloriaAll.csv');
lagloria(lagloria==0)=NaN;
for i=1:size(lagloria,2)
    lagloria(:,i) = lagloria(:,i) - min(lagloria(:,i));
    lagloria(:,i) = lagloria(:,i)./max(lagloria(:,i));
end

%%%%%%% [11] Golden Horn Batolith (Eddy et al. 2016)
GHB=csvread('GoldenHornAll.csv');
GHB(GHB==0)=NaN;
for i=1:size(GHB,2)
    GHB(:,i) = GHB(:,i) - min(GHB(:,i));
    GHB(:,i) = GHB(:,i)./max(GHB(:,i));
end

%%%%%%% [12] Turkey Creek Caldera and intracaldera intrusions (Deering et al. 2016)
% Turkey Creek Intracaldera intrusions
TurCrG = csvread('TurkeyCreekGranites.csv');
TurCrG(TurCrG==0)=NaN;
for i=1:size(TurCrG,2)
    TurCrG(:,i) = TurCrG(:,i) - min(TurCrG(:,i));
    TurCrG(:,i) = TurCrG(:,i)./max(TurCrG(:,i));
end
% Turkey Creek Caldera
TurCrV = csvread('TurkeyCreekVolcanics.csv');
TurCrV(TurCrV==0)=NaN;
for i=1:size(TurCrV,2)
    TurCrV(:,i) = TurCrV(:,i) - min(TurCrV(:,i));
    TurCrV(:,i) = TurCrV(:,i)./max(TurCrV(:,i));
end

%%%%%%% [13] Latir Volcanic field and intrusives (Tappa et al. 2013)
% Latir intrusives
LatG = csvread('LatirGranites.csv');
LatG(LatG==0)=NaN;
for i=1:size(LatG,2)
    LatG(:,i) = LatG(:,i) - min(LatG(:,i));
    LatG(:,i) = LatG(:,i)./max(LatG(:,i));
end
% Latir volcanics
LatV = csvread('LatirVolcanics.csv');
LatV(LatV==0)=NaN;
for i=1:size(LatV,2)
    LatV(:,i) = LatV(:,i) - min(LatV(:,i));
    LatV(:,i) = LatV(:,i)./max(LatV(:,i));
end

%%%%%%% [14] Kneeling Nun Tuff (Szymanowski et al., 2019)
KNT=csvread('KneelingNunAll.csv');
KNT(KNT==0)=NaN;
for i=1:size(KNT,2)
    KNT(:,i) = KNT(:,i) - min(KNT(:,i));
    KNT(:,i) = KNT(:,i)./max(KNT(:,i));
end

%%%%%%% Yellowstone Caldera
% [15] Huckleberry Ridge Tuff (Singer et al., 2014; Wotzlaw et al., 2015)
HRT=csvread('HuckleberryRidgenoANT.csv');
HRT(HRT==0)=NaN;
for i=1:size(HRT,2)
    HRT(:,i) = HRT(:,i) - min(HRT(:,i));
    HRT(:,i) = HRT(:,i)./max(HRT(:,i));
end
% [16] Lava Creek Tuff (Matthews et al., 2015; Wotzlaw et al., 2015)
LCT=csvread('LavaCreekAll.csv');
LCT(LCT==0)=NaN;
for i=1:size(LCT,2)
    LCT(:,i) = LCT(:,i) - min(LCT(:,i));
    LCT(:,i) = LCT(:,i)./max(LCT(:,i));
end

%%%%%%% [17] Fish Canyon Tuff (Wotzlaw et al., 2013)
FC = csvread('FishCanyonAll.csv'); 
FC(FC==0)=NaN;
for i=1:size(FC,2)
    FC(:,i) = FC(:,i) - min(FC(:,i));
    FC(:,i) = FC(:,i)./max(FC(:,i));
end



%%%%%%% Kernel Density Estimate (KDE) of averaged zircon age spectra, 
%%%%%%% scaled between initiation and termination of zircon crystallization
% Calculate density distributions
npoints = 100;
[bergell_f,bergell_x,bergell_bw] = ksdensity(bergell(:));
[elba_f,elba_x,elba_bw] = ksdensity(elba(:));
[MtGivens_f,MtGivens_x,MtGivens_bw] = ksdensity(MtGivens(:));
[MtStuart_f,MtStuart_x,MtStuart_bw] = ksdensity(MtStuart(:));
[TPeak_f,TPeak_x,TPeak_bw] = ksdensity(TPeak(:));
[JMuir_f,JMuir_x,JMuir_bw] = ksdensity(JMuir(:));
[Tuolomne_f,Tuolomne_x,Tuolomne_bw] = ksdensity(Tuolomne(:));
[LDVC_f,LDVC_x,LDVC_bw] = ksdensity(LDVC(:));
[VFIC_f,VFIC_x,VFIC_bw] = ksdensity(VFIC(:));
[lagloria_f,lagloria_x,lagloria_bw] = ksdensity(lagloria(:));
[GHB_f,GHB_x,GHB_bw] = ksdensity(GHB(:));
[TurCrG_f,TurCrG_x,TurCrG_bw] = ksdensity(TurCrG(:));
[TurCrV_f,TurCrV_x,TurCrV_bw] = ksdensity(TurCrV(:));
[LatG_f,LatG_x,LatG_bw] = ksdensity(LatG(:));
[LatV_f,LatV_x,LatV_bw] = ksdensity(LatV(:));
[SesiaG_f,SesiaG_x,SesiaG_bw] = ksdensity(SesiaG(:));
[SesiaV_f,SesiaV_x,SesiaV_bw] = ksdensity(SesiaV(:));
[SesiaGNA_f,SesiaGNA_x,SesiaGNA_bw] = ksdensity(SesiaGNA(:));
[SesiaVNA_f,SesiaVNA_x,SesiaVNA_bw] = ksdensity(SesiaVNA(:));
[KNT_f,KNT_x,KNT_bw] = ksdensity(KNT(:));
[HRT_f,HRT_x,HRT_bw] = ksdensity(HRT(:));
[LCT_f,LCT_x,LCT_bw] = ksdensity(LCT(:));
[FC_f,FC_x,FC_bw] = ksdensity(FC(:));
% Obtain ksdensity values between 0-bandwith and 1+bandwidth
bdist = interp1(bergell_x,bergell_f,linspace(0-0.05,1+bergell_bw,npoints));
edist = interp1(elba_x,elba_f,linspace(0-0.05,1+elba_bw,npoints));
MtGdist = interp1(MtGivens_x,MtGivens_f,linspace(0-0.05,1+MtGivens_bw,npoints));
TPeakdist = interp1(TPeak_x,TPeak_f,linspace(0-0.05,1+TPeak_bw,npoints));
MtSdist = interp1(MtStuart_x,MtStuart_f,linspace(0-0.05,1+MtStuart_bw,npoints));
Tuolomnedist = interp1(Tuolomne_x,Tuolomne_f,linspace(0-0.05,1+Tuolomne_bw,npoints));
JMuirdist = interp1(JMuir_x,JMuir_f,linspace(0-0.05,1+JMuir_bw,npoints));
LDVCdist = interp1(LDVC_x,LDVC_f,linspace(0-0.05,1+LDVC_bw,npoints));
VFICdist = interp1(VFIC_x,VFIC_f,linspace(0-0.05,1+VFIC_bw,npoints));
GHBdist = interp1(GHB_x,GHB_f,linspace(0-0.05,1+GHB_bw,npoints));
lgdist = interp1(lagloria_x,lagloria_f,linspace(0-0.05,1+lagloria_bw,npoints));
TuCGdist = interp1(TurCrG_x,TurCrG_f,linspace(0-0.05,1+TurCrG_bw,npoints));
TuCVdist = interp1(TurCrV_x,TurCrV_f,linspace(0-0.05,1+TurCrV_bw,npoints));
LatGdist = interp1(LatG_x,LatG_f,linspace(0-0.05,1+LatG_bw,npoints));
LatVdist = interp1(LatV_x,LatV_f,linspace(0-0.05,1+LatV_bw,npoints));
SGdist = interp1(SesiaG_x,SesiaG_f,linspace(0-0.05,1+SesiaG_bw,npoints));
SVdist = interp1(SesiaV_x,SesiaV_f,linspace(0-0.05,1+SesiaV_bw,npoints));
SGNAdist = interp1(SesiaGNA_x,SesiaGNA_f,linspace(0-0.05,1+SesiaGNA_bw,npoints));
SVNAdist = interp1(SesiaVNA_x,SesiaVNA_f,linspace(0-0.05,1+SesiaVNA_bw,npoints));
Kdist = interp1(KNT_x,KNT_f,linspace(0-0.05,1+KNT_bw,npoints));
HRTdist = interp1(HRT_x,HRT_f,linspace(0-0.05,1+HRT_bw,npoints));
LCTdist = interp1(LCT_x,LCT_f,linspace(0-0.05,1+LCT_bw,npoints));
FCdist = interp1(FC_x,FC_f,linspace(0-0.05*FC_bw,1+FC_bw,npoints));
% Normalize each distribution
x = linspace(0,1,npoints);
bdist = bdist./trapz(x,bdist); Bergell_pluton=bdist;
edist = edist./trapz(x,edist); Capanne_pluton=edist;
MtGdist = MtGdist./trapz(x,MtGdist); Mount_Givens_pluton=MtGdist;
MtSdist = MtSdist./trapz(x,MtSdist); Mount_Stuart_pluton=MtSdist;
TPeakdist = TPeakdist./trapz(x,TPeakdist); Tenpeak_pluton=TPeakdist;
JMuirdist = JMuirdist./trapz(x,JMuirdist); JohnMuir_IntrusiveSuite=JMuirdist;
Tuolomnedist = Tuolomnedist./trapz(x,Tuolomnedist); Tuolomne_IntrusiveSuite=Tuolomnedist;
LDVCdist = LDVCdist./trapz(x,LDVCdist); LagodellaVacca_IntrusiveComplex=LDVCdist;
VFICdist = VFICdist./trapz(x,VFICdist); ValFredda_IntrusiveComplex=VFICdist;
lgdist = lgdist./trapz(x,lgdist); LaGloria_pluton=lgdist;
GHBdist = GHBdist./trapz(x,GHBdist); Golden_Horn_Batholith=GHBdist;
TuCGdist = TuCGdist./trapz(x,TuCGdist); Turkey_Creek_Intrusives=TuCGdist;
TuCVdist = TuCVdist./trapz(x,TuCVdist); Turkey_Creek_Caldera=TuCVdist;
LatGdist = LatGdist./trapz(x,LatGdist); Latir_Intrusives=LatGdist;
LatVdist = LatVdist./trapz(x,LatVdist); Latir_Volcanics=LatVdist;
SGdist = SGdist./trapz(x,SGdist); Valle_Mosso_Pluton=SGdist;
SVdist = SVdist./trapz(x,SVdist); Sesia_Caldera=SVdist;
SGNAdist = SGNAdist./trapz(x,SGNAdist); Valle_Mosso_Pluton_NoAntecrysts=SGNAdist;
SVNAdist = SVNAdist./trapz(x,SVNAdist); Sesia_Caldera_NoAntecrysts=SVNAdist;
Kdist = Kdist./trapz(x,Kdist); Kneeling_Nun_Tuff=Kdist;
HRTdist = HRTdist./trapz(x,HRTdist); Huckleberry_Ridge_Tuff=HRTdist;
LCTdist = LCTdist./trapz(x,LCTdist); Lava_Creek_Tuff=LCTdist;
FCdist = FCdist./trapz(x,FCdist); Fish_Canyon_Tuff=FCdist;


                    %%%%%% [User Input] %%%%%% 
%%% Load from .csv file and calculate KDE for new zircon ages distribution 
%%% by editing the mock script:

% ExampleName = csvread('Example.csv'); 
% ExampleName(ExampleName==0)=NaN;
% for i=1:size(ExampleName,2)
%    ExampleName(:,i) = ExampleName(:,i) - min(ExampleName(:,i));
%    ExampleName(:,i) = ExampleName(:,i)./max(ExampleName(:,i));
% end
% [EX_f,EX_x,EX_bw] = ksdensity(ExampleName(:));
% EXdist = interp1(EX_x,EX_f,linspace(0-0.05*EX_bw,1+EX_bw,npoints));
% EXdist = EXdist./trapz(x,EXdist); Example_Name=EXdist;


% For questions and ideas on how to improve the code please contact: 
% ltavazzani@smu.edu
