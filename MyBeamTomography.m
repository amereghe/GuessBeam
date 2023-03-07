% {}~

%% include libraries
% - include Matlab libraries
pathToLibrary="../MatLabTools";
addpath(genpath(pathToLibrary));
% - include lib libraries
pathToLibrary="lib";
addpath(genpath(pathToLibrary));

%%
% parse summary file and get phase advance

whatScan="MUX";
eleScan="X3_011B_VWN";
iFileName=sprintf("../optics/HEBT/output/summary_%s_%s.tfs",eleScan,whatScan);
myTitle=sprintf("%s scan at %s - ISO3 (XPR) - Proton - 30mm",whatScan,LabelMe(eleScan));

summdaryData=readmatrix(iFileName,'HeaderLines',1,'Delimiter',',','FileType','text');
if (strcmpi(whatScan,"MUX"))
    MUs=summdaryData(:,end-1);
elseif (strcmpi(whatScan,"MUY"))
    MUs=summdaryData(:,end);
else
    clear MUs;
end
MUs=rad2deg(MUs*2*pi);

%% define bar of charge

% bar of charge dimensions
bb=5E-3;    % [m]
hh=0.1E-3;  % [rad]

% sample points
nPoints=10000;
MyPointsStart=BeamSample_Rect(nPoints,bb,hh);

% contours
clear MyContoursStart; MyContoursStart=missing();
% contouring Rect
MyContoursStart=ExpandMat(MyContoursStart,GenPointsAlongRectangle(bb,hh));
% contouring ellypse
clear alphaStart betaStart emiGstart;
[alphaStart,betaStart,emiGstart]=GetOpticsFromSigmaMatrix(MyPointsStart);   % ellypse orientation
emiGstart=max(GetSinglePartEmittance(MyPointsStart,alphaStart,betaStart));  % max 
MyContoursStart=ExpandMat(MyContoursStart,GenPointsAlongEllypse(alphaStart,betaStart,emiGstart));

% show starting configuration
% % - just points and contours
% ShowMe(MyPointsStart,MyContoursStart);
% - 2D histograms
clear nCountsStart nCountsXstart nCountsYstart xbStart ybStart
[nCountsStart,nCountsXstart,nCountsYstart,xbStart,ybStart]=Get2dHistograms(MyPointsStart(:,1),MyPointsStart(:,2));
% figure();
% Plot2DHistograms(nCountsStart,nCountsXstart,nCountsYstart,xbStart,ybStart,"z [m]","zp [rad]",MyContoursStart);
% sgtitle("starting bar of charge");
% - points and 1D histograms
figure();
Plot2DHistograms(MyPointsStart,nCountsXstart,nCountsYstart,xbStart,ybStart,"z [m]","zp [rad]",MyContoursStart,false);
sgtitle("starting bar of charge");

%% rotate distribution
angle=MUs(1);
RM=Rot2D(-angle);
% % - in normalised phase space
% lNorm=true;
% MyPointsStartTransported=AdvanceMyPoints(MyPointsStart,RM,lNorm);
% - in physical phase space
lNorm=false;
MyPointsStartTransported=AdvanceMyPoints(MyPointsStart,RM,lNorm);

%% rotate Bar

lNorm=false;
angles=MUs(2:end);
clear RMs; RMs=Rot2D(-angles);
rr=max(bb,hh)*1.5; % [mm]
xb=linspace(-rr,rr,100+1)*1E3; % [m] to [mm]
yb=linspace(-rr,rr,100+1)*1E3; % [rad] to [mrad]

clear alphaFit betaFit emiGfit phisFit; alphaFit=NaN(size(RMs,3),1); betaFit=NaN(size(RMs,3),1); emiGfit=NaN(size(RMs,3),1); phisFit=NaN(size(RMs,3),1);
clear alphaPrp betaPrp emiGprp phisPrp; alphaPrp=NaN(size(RMs,3),1); betaPrp=NaN(size(RMs,3),1); emiGprp=NaN(size(RMs,3),1); phisPrp=NaN(size(RMs,3),1);
clear sigZ sigZP; sigZ=NaN(size(RMs,3),1); sigZP=NaN(size(RMs,3),1);
clear profiles; profiles=NaN(length(xb)-1,1+size(RMs,3));
[alphaFit,betaFit,emiGfit,phisFit,alphaPrp,betaPrp,emiGprp,phisPrp,sigZ,sigZP,profiles]=...
    ScanBar(RMs,MyPointsStart,xb,yb,betaStart,alphaStart,lNorm,"Rotating Bar");
Show3Ddistributions(profiles,angles,"\theta [deg]","","Rotating Bar");
% ShowScans(angles,"\theta [deg]",alphaFit,betaFit,emiGfit,phisFit,sigZ,sigZP,"","Rotating Bar");
% ComparePropagationVsFit(angles,"\theta [deg]",alphaFit,betaFit,emiGfit,phisFit,alphaPrp,betaPrp,emiGprp,phisPrp,sigZ,sigZP,"Rotating Bar");
