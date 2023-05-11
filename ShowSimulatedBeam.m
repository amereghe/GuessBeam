% {}~

%% include libraries
% - include Matlab libraries
pathToLibrary="../MatLabTools";
addpath(genpath(pathToLibrary));
% - include lib libraries
pathToLibrary="lib";
addpath(genpath(pathToLibrary));

%% parse optics file
whatScan="MUX";
betScan=sprintf("BET%s",extractBetween(whatScan,strlength(whatScan),strlength(whatScan)));
alfScan=sprintf("ALF%s",extractBetween(whatScan,strlength(whatScan),strlength(whatScan)));
% starting optics functions: V.Lante, "XPR optics design"
betaStartNom=8.184; alfaStartNom=-0.443; gammaStartNom=(1+alfaStartNom^2)/betaStartNom;
emiGeo=5E-6; % [m rad]
sigZStartNom=sqrt(betaStartNom*emiGeo); sigZPStartNom=sqrt(gammaStartNom*emiGeo);
iFileNameOpticsSummary="../optics/HEBT/output_p030_MUX_X3-011B-VWN_free_10p0degs/x3_011b_vwn_summary_optics.tfs";
iFileNameRMatrixSummary="../optics/HEBT/output_p030_MUX_X3-011B-VWN_free_10p0degs/x3_011b_vwn_summary_rmatrix.tfs";
% - crunch optics file
scannedOptics=ParseTfsTable(iFileNameOpticsSummary,'optics',"SCAN");
[ colNames, colUnits, colFacts, mapping, readFormat ] = ...
        GetColumnsAndMappingTFS("optics","SCAN");
iID=mapping(find(strcmp(colNames,'ID')));
iMU=mapping(find(strcmp(colNames,whatScan)));
iMUX=mapping(find(strcmp(colNames,"MUX")));
iMUY=mapping(find(strcmp(colNames,"MUY")));
zeroPos=find([scannedOptics{:,iID}' 0]==0);
sortedIDs=NaN(size(size(scannedOptics{1}),1),1); iLast=0;
for ii=1:length(zeroPos)-1
    myRange=zeroPos(ii):zeroPos(ii+1)-1;
    [~,mySortedIDs]=sort(scannedOptics{iMU}(myRange)); % sort by phase advance
    nScanPoints=length(mySortedIDs);
    sortedIDs(1+iLast:nScanPoints+iLast)=mySortedIDs+iLast;
    iLast=iLast+nScanPoints;
end
nScanPoints=length(sortedIDs);
% - get reference optics
mu0=rad2deg(scannedOptics{iMU}(zeroPos(1))*2*pi);
beta0=scannedOptics{mapping(find(strcmp(colNames,betScan)))}(zeroPos(1));
alfa0=scannedOptics{mapping(find(strcmp(colNames,alfScan)))}(zeroPos(1));
gamma0=(1+alfa0^2)/beta0;
% - crunch response matrix file
scannedRMatrices=ParseTfsTable(iFileNameRMatrixSummary,'RMatrix',"SCAN");
[ colNames, colUnits, colFacts, mapping, readFormat ] = ...
        GetColumnsAndMappingTFS("RMatrix","SCAN");
readRMs=GetTransMatrix(cell2mat(scannedRMatrices),"HOR","SCAN");
MUs=rad2deg(scannedOptics{iMU}*2*pi);
[uMUs,iuMUs,~]=unique(MUs(sortedIDs));
nActualPoints=length(MUs(sortedIDs(iuMUs)));

%% define bar of charge

% bar of charge dimensions
bb=5E-3;    % [m]
hh=0.1E-3;  % [rad]

% sample points
nParticles=10000;
MyPointsStart=BeamSample_Rect(nParticles,bb,hh);

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
% % - 2D histograms (physical units)
% clear nCountsStart nCountsXstart nCountsYstart xbStart ybStart
% [nCountsStart,nCountsXstart,nCountsYstart,xbStart,ybStart]=Get2dHistograms(MyPointsStart(:,1),MyPointsStart(:,2));
% % - 2D histograms (normalised units)
% clear nCountsStartNorm nCountsXstartNorm nCountsYstartNorm xbStartNorm ybStartNorm
% MyPointsStartNorm=Phys2Norm(MyPointsStart,betaStartNom,alfaStartNom,emiGeo);
% [nCountsStartNorm,nCountsXstartNorm,nCountsYstartNorm,xbStartNorm,ybStartNorm]=Get2dHistograms(MyPointsStartNorm(:,1),MyPointsStartNorm(:,2));
% % figure();
% % Plot2DHistograms(nCountsStart,nCountsXstart,nCountsYstart,xbStart,ybStart,"z [m]","zp [rad]",MyContoursStart);
% % sgtitle("starting bar of charge");
% % - points and 1D histograms
% figure();
% Plot2DHistograms(MyPointsStart,nCountsXstart,nCountsYstart,xbStart,ybStart,"z [m]","zp [rad]",MyContoursStart,false);
% sgtitle("starting bar of charge");
% 
%% show 2D histogram
ShowPhysNorm(MyPointsStart,betaStartNom,alfaStartNom,emiGeo,"starting beam");

%% transport bar
lNorm=false;
angles=MUs(sortedIDs(iuMUs));
% clear RMs; RMs=Rot2D(-angles);
RMs=readRMs(1:2,1:2,sortedIDs(iuMUs));
rx=50E-03; % [m]
ry=3E-03; % [rad]
xb=linspace(-rx,rx,100+1)*1E3; % [m] to [mm]
yb=linspace(-ry,ry,100+1)*1E3; % [rad] to [mrad]

clear alphaFit betaFit emiGfit phisFit; alphaFit=NaN(size(RMs,3),1); betaFit=NaN(size(RMs,3),1); emiGfit=NaN(size(RMs,3),1); phisFit=NaN(size(RMs,3),1);
clear alphaPrp betaPrp emiGprp phisPrp; alphaPrp=NaN(size(RMs,3),1); betaPrp=NaN(size(RMs,3),1); emiGprp=NaN(size(RMs,3),1); phisPrp=NaN(size(RMs,3),1);
clear sigZ sigZP; sigZ=NaN(size(RMs,3),1); sigZP=NaN(size(RMs,3),1);
clear profiles; profiles=NaN(length(xb)-1,1+size(RMs,3));
[alphaFit,betaFit,emiGfit,phisFit,alphaPrp,betaPrp,emiGprp,phisPrp,sigZ,sigZP,profiles]=...
    ScanBar(RMs,MyPointsStart,xb,yb,betaStart,alphaStart,lNorm,"Rotating Bar");
Show3Ddistributions(profiles,angles,"\theta [deg]","","Rotating Bar");
% ShowScans(angles,"\theta [deg]",alphaFit,betaFit,emiGfit,phisFit,sigZ,sigZP,"","Rotating Bar");
% ComparePropagationVsFit(angles,"\theta [deg]",alphaFit,betaFit,emiGfit,phisFit,alphaPrp,betaPrp,emiGprp,phisPrp,sigZ,sigZP,"Rotating Bar");

% try to get inverse radon transform
ShowSinogramIR(profiles,angles,1:nActualPoints,sqrt(gamma0/beta0));

%% transport beam through nominal optics
RM=RMs(:,:,mu0==angles);
% - in physical phase space (lNorm=true for normalised units)
lNorm=false;
MyPointsTransportedNominal=AdvanceMyPoints(MyPointsStart,RM,lNorm);

% show 2D histogram
ShowPhysNorm(MyPointsTransportedNominal,beta0,alfa0,emiGeo,"nominal transported beam");


%% scan angular interval
iMin=0;
for iMax=find(diff(angles>=angles(1)+189)):length(angles)
    iMin=iMin+1;
    myRange=iMin:iMax;
    IR = iradon(profiles(:,1+myRange),angles(myRange));
    % figure();
    imshow(IR,'InitialMagnification',400); colormap(hot); 
    title(sprintf("\\theta[degs]=[%g:%g] - \\Delta\\theta[degs]=%g",angles(iMin),angles(iMax),angles(iMax)-angles(iMin)));
    pause(1);
    % close gcf;
end

%% local functions
function ShowPhysNorm(MyPoints,beta,alfa,emiGeo,myTitle,lPhys)

    if (~exist("lPhys","var")), lPhys=true; end
    if (lPhys)
        MyPointsPhys=MyPoints;
        MyPointsNorm=Phys2Norm(MyPointsPhys,beta,alfa,emiGeo);
    else
        MyPointsNorm=MyPoints;
        MyPointsPhys=Norm2Phys(MyPointsNorm,beta,alfa,emiGeo);
    end
    
    %% compute histograms
    clear nCountsPhys nCountsXphys nCountsYphys xbPhys ybPhys;
    [nCountsPhys,nCountsXphys,nCountsYphys,xbPhys,ybPhys]=Get2dHistograms(MyPointsPhys(:,1),MyPointsPhys(:,2));
    clear nCountsNorm nCountsXnorm nCountsYnorm xbNorm ybNorm;
    [nCountsNorm,nCountsXnorm,nCountsYnorm,xbNorm,ybNorm]=Get2dHistograms(MyPointsNorm(:,1),MyPointsNorm(:,2));
    
    %% get ready for plotting
    % - physical units
    XsPhys=xbPhys(2:end)-diff(xbPhys(1:2))/2;
    YsPhys=ybPhys(2:end)-diff(ybPhys(1:2))/2;
    % - normalised units
    XsNorm=xbNorm(2:end)-diff(xbNorm(1:2))/2;
    YsNorm=ybNorm(2:end)-diff(ybNorm(1:2))/2;
    
    %% actually show
    figure();
    % - physical units
    subplot(1,2,1);
    imagesc('XData',XsPhys,'YData',YsPhys,'CData',nCountsPhys);
    grid(); colorbar;
    xlabel("z [m]"); ylabel("zp [rad]"); title("physical units");
    % - normalised units
    subplot(1,2,2);
    imagesc('XData',XsNorm,'YData',YsNorm,'CData',nCountsNorm);
    grid(); colorbar;
    xlabel("z []"); ylabel("zp []"); title("normalised units");
    %
    sgtitle(myTitle);
    
end

function ShowSinogramIR(profiles,angles,myRange,Ryx)
    if (~exist("myRange","var")), myRange=1:length(angles); end
    if (~exist("Ryx","var")), Ryx=1; end

    %% compute Radon anti-transform
    IR = iradon(profiles(:,1+myRange),angles(myRange),"spline","Hamming");

    %% actually plot  
    figure();
    % - sinogram
    subplot(1,2,1);
    imagesc('XData',angles,'YData',profiles(:,1),'CData',profiles(:,2:end));
    grid(); colorbar;
    xlabel('\mu [degs]'); ylabel('z [mm]'); title('Synogram');
    % - inverse radon transform
    xs=linspace(min(profiles(:,1)),max(profiles(:,1)),size(IR,1)); % [mm]
    ys=xs*Ryx; % [mrad]
    subplot(1,2,2);
    imagesc('XData',xs,'YData',ys,'CData',IR);
    grid(); colorbar;
    xlabel('z [mm]'); ylabel('zp [mrad]'); title('reconstructed');
end
