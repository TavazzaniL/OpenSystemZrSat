% Data from literature used in bootstrapped crystallization dist. plots
cd /Users/Utente/Desktop/zurich_ID-TIMS/Keller_papers_zircSat/BayeZirChron.c-master/examples/BootstrappedCrystallizationDistributionExamples/UPB_LIT

% Bergell pluton (Samperton et al.2015)
bergell=csvread('BergellAll.csv');
bergell(bergell==0)=NaN;
for i=1:size(bergell,2)
    bergell(:,i) = bergell(:,i) - min(bergell(:,i));
    bergell(:,i) = bergell(:,i)./max(bergell(:,i));
end

% Mt Capanne pluton (Barboni et al.2015)
elba = csvread('ElbaAll.csv');
elba(elba==0)=NaN;
for i=1:size(elba,2)
    elba(:,i) = elba(:,i) - min(elba(:,i));
    elba(:,i) = elba(:,i)./max(elba(:,i));
end

%Mt Givens granodiorite (Frazer et al., 2015)
MtGivens = csvread('MTGIVENS.csv');
MtGivens(MtGivens==0)=NaN;
for i=1:size(MtGivens,2)
    MtGivens(:,i) = MtGivens(:,i) - min(MtGivens(:,i));
    MtGivens(:,i) = MtGivens(:,i)./max(MtGivens(:,i));
end

%Mt Stuart Intrusive complex and Tenpeak pluton (Matzel et al., 2007)
MtStuart = csvread('MTSTUART.csv');
MtStuart(MtStuart==0)=NaN;
for i=1:size(MtStuart,2)
    MtStuart(:,i) = MtStuart(:,i) - min(MtStuart(:,i));
    MtStuart(:,i) = MtStuart(:,i)./max(MtStuart(:,i));
end
TPeak = csvread('Tenpeak.csv');
TPeak(TPeak==0)=NaN;
for i=1:size(TPeak,2)
    TPeak(:,i) = TPeak(:,i) - min(TPeak(:,i));
    TPeak(:,i) = TPeak(:,i)./max(TPeak(:,i));
end

%Tuolomne Intrusive suite (Coleman et al., 2004)
Tuolomne = csvread('Tuolomne.csv');
Tuolomne(Tuolomne==0)=NaN;
for i=1:size(Tuolomne,2)
    Tuolomne(:,i) = Tuolomne(:,i) - min(Tuolomne(:,i));
    Tuolomne(:,i) = Tuolomne(:,i)./max(Tuolomne(:,i));
end

%John Muir Intrusive suite (Davis et al., 2011)
JMuir = csvread('JMuir.csv');
JMuir(JMuir==0)=NaN;
for i=1:size(JMuir,2)
    JMuir(:,i) = JMuir(:,i) - min(JMuir(:,i));
    JMuir(:,i) = JMuir(:,i)./max(JMuir(:,i));
end

%Southern Adamello (Lago della Vacca, Schoene et al., 2012; Val Fredda, Broderick et al., 2015)
LDVC = csvread('LDVC.csv');
LDVC(LDVC==0)=NaN;
for i=1:size(LDVC,2)
    LDVC(:,i) = LDVC(:,i) - min(LDVC(:,i));
    LDVC(:,i) = LDVC(:,i)./max(LDVC(:,i));
end
VFIC = csvread('VFIC.csv');
VFIC(VFIC==0)=NaN;
for i=1:size(VFIC,2)
    VFIC(:,i) = VFIC(:,i) - min(VFIC(:,i));
    VFIC(:,i) = VFIC(:,i)./max(VFIC(:,i));
end

%La Gloria pluton (Gutierrez et al. 2018)
lagloria=csvread('LAGLORIA.csv');
lagloria(lagloria==0)=NaN;
for i=1:size(lagloria,2)
    lagloria(:,i) = lagloria(:,i) - min(lagloria(:,i));
    lagloria(:,i) = lagloria(:,i)./max(lagloria(:,i));
end

%Golden Horn Batolith (Eddy et al. 2018)
GHB=csvread('GHB.csv');
GHB(GHB==0)=NaN;
for i=1:size(GHB,2)
    GHB(:,i) = GHB(:,i) - min(GHB(:,i));
    GHB(:,i) = GHB(:,i)./max(GHB(:,i));
end

%Turkey Creek Caldera and intracaldera intrusions (Deering et al. 2016)
TurCrG = csvread('TurCrG.csv');
TurCrG(TurCrG==0)=NaN;
for i=1:size(TurCrG,2)
    TurCrG(:,i) = TurCrG(:,i) - min(TurCrG(:,i));
    TurCrG(:,i) = TurCrG(:,i)./max(TurCrG(:,i));
end
TurCrV = csvread('TurCrV.csv');
TurCrV(TurCrV==0)=NaN;
for i=1:size(TurCrV,2)
    TurCrV(:,i) = TurCrV(:,i) - min(TurCrV(:,i));
    TurCrV(:,i) = TurCrV(:,i)./max(TurCrV(:,i));
end

%Latir Volcanic field and intrusives (Tappa et al. 2013)
LatG = csvread('LatG.csv');
LatG(LatG==0)=NaN;
for i=1:size(LatG,2)
    LatG(:,i) = LatG(:,i) - min(LatG(:,i));
    LatG(:,i) = LatG(:,i)./max(LatG(:,i));
end
LatV = csvread('LatV.csv');
LatV(LatV==0)=NaN;
for i=1:size(LatV,2)
    LatV(:,i) = LatV(:,i) - min(LatV(:,i));
    LatV(:,i) = LatV(:,i)./max(LatV(:,i));
end

SesiaG=csvread('SesiaG1.csv');
SesiaG(SesiaG==0)=NaN;
for i=1:size(SesiaG,2)
    SesiaG(:,i) = SesiaG(:,i) - min(SesiaG(:,i));
    SesiaG(:,i) = SesiaG(:,i)./max(SesiaG(:,i));
end
SesiaV=csvread('SesiaV.csv');
SesiaV(SesiaV==0)=NaN;
for i=1:size(SesiaV,2)
    SesiaV(:,i) = SesiaV(:,i) - min(SesiaV(:,i));
    SesiaV(:,i) = SesiaV(:,i)./max(SesiaV(:,i));
end

SesiaGNA=csvread('SesiaG_NoANT1.csv');
SesiaGNA(SesiaGNA==0)=NaN;
for i=1:size(SesiaGNA,2)
    SesiaGNA(:,i) = SesiaGNA(:,i) - min(SesiaGNA(:,i));
    SesiaGNA(:,i) = SesiaGNA(:,i)./max(SesiaGNA(:,i));
end
SesiaVNA=csvread('SesiaV_NoANT.csv');
SesiaVNA(SesiaVNA==0)=NaN;
for i=1:size(SesiaVNA,2)
    SesiaVNA(:,i) = SesiaVNA(:,i) - min(SesiaVNA(:,i));
    SesiaVNA(:,i) = SesiaVNA(:,i)./max(SesiaVNA(:,i));
end

%Kneeling Nun Tuff (Szymanowski et al., 2019)
KNT=csvread('KNTAll.csv');
KNT(KNT==0)=NaN;
for i=1:size(KNT,2)
    KNT(:,i) = KNT(:,i) - min(KNT(:,i));
    KNT(:,i) = KNT(:,i)./max(KNT(:,i));
end

%Hucklebarry Ridge Tuff (Singer et al., 2011; Wotzlaw et al., 2015)
HRT=csvread('HRT_noANT.csv');
HRT(HRT==0)=NaN;
for i=1:size(HRT,2)
    HRT(:,i) = HRT(:,i) - min(HRT(:,i));
    HRT(:,i) = HRT(:,i)./max(HRT(:,i));
end

LCT=csvread('LCT.csv');
LCT(LCT==0)=NaN;
for i=1:size(LCT,2)
    LCT(:,i) = LCT(:,i) - min(LCT(:,i));
    LCT(:,i) = LCT(:,i)./max(LCT(:,i));
end

%Fish Canyon Tuff (Wotzlaw et al., 2013)
FC = csvread('FishCanyonAll.csv'); 
FC(FC==0)=NaN;
for i=1:size(FC,2)
    FC(:,i) = FC(:,i) - min(FC(:,i));
    FC(:,i) = FC(:,i)./max(FC(:,i));
end

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

% Obtain ksdensity values between 0-bandwith and 1+2*bandwidth
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
bdist = bdist./trapz(x,bdist);
edist = edist./trapz(x,edist);
MtGdist = MtGdist./trapz(x,MtGdist);
MtSdist = MtSdist./trapz(x,MtSdist);
TPeakdist = TPeakdist./trapz(x,TPeakdist);
JMuirdist = JMuirdist./trapz(x,JMuirdist);
Tuolomnedist = Tuolomnedist./trapz(x,Tuolomnedist);
LDVCdist = LDVCdist./trapz(x,LDVCdist);
VFICdist = VFICdist./trapz(x,VFICdist);
lgdist = lgdist./trapz(x,lgdist);
GHBdist = GHBdist./trapz(x,GHBdist);
TuCGdist = TuCGdist./trapz(x,TuCGdist);
TuCVdist = TuCVdist./trapz(x,TuCVdist);
LatGdist = LatGdist./trapz(x,LatGdist);
LatVdist = LatVdist./trapz(x,LatVdist);
SGdist = SGdist./trapz(x,SGdist);
SVdist = SVdist./trapz(x,SVdist);
SGNAdist = SGNAdist./trapz(x,SGNAdist);
SVNAdist = SVNAdist./trapz(x,SVNAdist);
Kdist = Kdist./trapz(x,Kdist);
HRTdist = HRTdist./trapz(x,HRTdist);
LCTdist = LCTdist./trapz(x,LCTdist);
FCdist = FCdist./trapz(x,FCdist);

