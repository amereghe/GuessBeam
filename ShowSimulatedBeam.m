% {}~

%% include libraries
% % - include Matlab libraries
% pathToLibrary="../MatLabTools";
% addpath(genpath(pathToLibrary));
% % - include lib libraries
% pathToLibrary="lib";
% addpath(genpath(pathToLibrary));

%% user settings
whatScan="MUX";

% - starting optics functions: V.Lante, "XPR optics design"
BETXstartNom=8.184; ALFXstartNom=-0.443;
BETYstartNom=4.526; ALFYstartNom=-1.985;
emiGeo=1; % [m rad]

% - path to files
% iFileNameOpticsSummary="../optics/HEBT/output_p030_MUX_X3-011B-VWN_free_10p0degs/x3_011b_vwn_summary_optics.tfs";
% iFileNameRMatrixSummary="../optics/HEBT/output_p030_MUX_X3-011B-VWN_free_10p0degs/x3_011b_vwn_summary_rmatrix.tfs";
% iFileNameOpticsSummary="../optics/HEBT/output_p030_MUY_X3-011B-VWN_free_10p0degs/x3_011b_vwn_summary_optics.tfs";
% iFileNameRMatrixSummary="../optics/HEBT/output_p030_MUY_X3-011B-VWN_free_10p0degs/x3_011b_vwn_summary_rmatrix.tfs";
% iFileNameOpticsSummary="../optics/HEBT/output_p030_MUX_X3-011B-VWN_constrained_10p0degs/x3_011b_vwn_summary_optics.tfs";
% iFileNameRMatrixSummary="../optics/HEBT/output_p030_MUX_X3-011B-VWN_constrained_10p0degs/x3_011b_vwn_summary_rmatrix.tfs";
iFileNameOpticsSummary="../optics/HEBT/output_p030_MUX_X3-011B-VWN_constrained_02p5degs/x3_011b_vwn_summary_optics.tfs";
iFileNameRMatrixSummary="../optics/HEBT/output_p030_MUX_X3-011B-VWN_constrained_02p5degs/x3_011b_vwn_summary_rmatrix.tfs";

% - 2D histograms (Poincaré phase plots)
rx=64E-03; % used also for cameretta [m]
ry=3E-03;  % [rad]

% - beam sampling
nParticles=10000;

% - bar of charge dimensions
bb=20E-3;    % [sqrt(m)] or [m]
hh=0.1E-3;  % [sqrt(m)] or [rad]
barIsNormal=false; % bar of charge sampled in normal coordinates

% - Gaussian sampling
GaussIsNormal=true;

%% do the job

% - prepare variables
planeObs=extractBetween(whatScan,strlength(whatScan),strlength(whatScan));
betScan=sprintf("BET%s",planeObs);
alfScan=sprintf("ALF%s",planeObs);
if (strcmpi(planeObs,"X"))
    betaStartNom=BETXstartNom; alfaStartNom=ALFXstartNom;
    lBar=true; % track bar of charge
    scanPlane="HOR";
elseif (strcmpi(planeObs,"Y"))
    betaStartNom=BETYstartNom; alfaStartNom=ALFYstartNom;
    lBar=false; % track Gaussian beam
    scanPlane="VER";
end
gammaStartNom=(1+alfaStartNom^2)/betaStartNom;

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
readRMs=GetTransMatrix(cell2mat(scannedRMatrices),scanPlane,"SCAN");
MUs=rad2deg(scannedOptics{iMU}*2*pi);
[uMUs,iuMUs,~]=unique(MUs(sortedIDs));
nActualPoints=length(MUs(sortedIDs(iuMUs)));

%% define beam

% sample points
if (lBar)
    MyPointsStart=BeamSample_Rect(nParticles,bb,hh);
    if (barIsNormal)
        MyPointsStart=Norm2Phys(MyPointsStart,betaStartNom,alfaStartNom,emiGeo);
    end
else
    MyPointsStart=BeamSample_Gauss(nParticles,alfaStartNom,betaStartNom,emiGeo);
    if (GaussIsNormal)
        MyPointsStart=MyPointsStart/1000;
    end
end

% contours
clear MyContoursStart; MyContoursStart=missing();
if (lBar)
    MyContoursStart=ExpandMat(MyContoursStart,GenPointsAlongRectangle(bb,hh));
else
    MyContoursStart=GenPointsAlongEllypse(alfaStartNom,betaStartNom,emiGeo);
end

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
% sgtitle("starting beam");

% show 2D histogram
ShowCoordsPhysNorm(MyPointsStart,betaStartNom,alfaStartNom,emiGeo,"starting beam");

%% scan beam,
lNorm=false;
angles=MUs(sortedIDs(iuMUs));
% clear RMs; RMs=Rot2D(-angles);
RMs=readRMs(1:2,1:2,sortedIDs(iuMUs));
xb=linspace(-rx,rx,128+1)*1E3; % [m] to [mm]
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

% show simulated profiles as a sinogram
ShowSinogramProfiles(profiles,angles);

%% transport beam through nominal optics
RM=RMs(:,:,mu0==angles);
lNorm=false; % lNorm=true for normalised units
MyPointsTransportedNominal=AdvanceMyPoints(MyPointsStart,RM,lNorm);
% show 2D histogram
ShowCoordsPhysNorm(MyPointsTransportedNominal,beta0,alfa0,emiGeo,"nominal transported beam");

%% try to get inverse radon transform
% - select range of points
myRange=1:nActualPoints;
% - manipulate profiles:
%   selected range
myProfiles=profiles(:,[1 myRange+1]);
%   normalised coordinates
myProfiles(:,1)=myProfiles(:,1)/(sqrt(beta0*emiGeo)*1E3);
% - manipulate angles:
%   mu=-Radon_angle
%   subtract mu0 if you want to compare to transported beam
myAngles=-(angles(myRange));
ShowSinogramIR(myProfiles,myAngles);

%% stop
return

%% scan angular interval
iMin=0; Delta=180;
for iMax=find(diff(angles>angles(1)+Delta)):length(angles)
    iMin=iMin+1;
    myRange=iMin:iMax;
    myTitle=sprintf("\\theta[degs]=[%g:%g] - \\Delta\\theta[degs]=%g",angles(iMin),angles(iMax),angles(iMax)-angles(iMin));
    ShowSinogramIR(myProfiles(:,[1 myRange+1]),myAngles(myRange),myTitle);
    pause(1);
end

%% local functions
function ShowCoordsPhysNorm(MyPoints,beta,alfa,emiGeo,myTitle,lPhys,l3D)

    if (~exist("l3D","var")), l3D=true; end
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
    ff=figure();
    ff.Position(3)=2*ff.Position(3);
    % - physical units
    subplot(1,2,1);
    if (l3D)
        % 3D plot
        surf(XsPhys,YsPhys,nCountsPhys',"EdgeColor","none");
        zlabel("counts []");
        view(2); % starts appearing as a 2D plot
    else
        % 2D plot
        imagesc('XData',XsPhys,'YData',YsPhys,'CData',nCountsPhys');
        grid(); colorbar;
    end
    xlabel("z [m]"); ylabel("zp [rad]"); title("physical units");
    % - normalised units
    subplot(1,2,2);
    if (l3D)
        % 3D plot
        surf(XsNorm,YsNorm,nCountsNorm',"EdgeColor","none");
        zlabel("counts []");
        view(2); % starts appearing as a 2D plot
    else
        % 2D plot
        imagesc('XData',XsNorm,'YData',YsNorm,'CData',nCountsNorm');
        grid(); colorbar;
    end
    zMin=min([xlim ylim]); zMax=max([xlim ylim]);
    xlim([zMin zMax]); ylim([zMin zMax]);
    axis square;
    if (emiGeo==1)
        xlabel("z [\surd{m}]"); ylabel("zp [\surd{m}]");
    else
        xlabel("z []"); ylabel("zp []");
    end
    title("normalised units");
    %
    sgtitle(myTitle);
    
end

function ShowSinogramIR(profiles,angles,myTitle,l3D)

    if (~exist("l3D","var")), l3D=true; end

    %% compute Radon anti-transform
    myInterpolation="spline"; myFilter="Hamming"; % "Hamming","Ram-Lak"
    IR = iradon(profiles(:,2:end),angles,myInterpolation,myFilter);

    %% actually show
    ff=figure();
    ff.Position(3)=2*ff.Position(3);
    % - sinogram
    subplot(1,2,1);
    PlotSinogramProfiles(profiles,angles);
    grid(); colorbar;
    xlabel('Angle [degs]'); ylabel("z [\surd{m}]"); title('sinogram');
    % - inverse radon transform (normalised coordinates)
    rs=linspace(min(profiles(:,1))/sqrt(2),max(profiles(:,1))/sqrt(2),size(IR,1)); % []
    subplot(1,2,2);
    if (l3D)
        % 3D plot
        surf(rs,rs,IR,"EdgeColor","none");
        zlabel("inverse Radon transform");
        view(2); % starts appearing as a 2D plot
    else
        % 2D plot
        imagesc('XData',rs,'YData',rs,'CData',IR);
        grid(); colorbar;
    end
    axis square;
    xlabel('z [\surd{m}]'); ylabel('zp [\surd{m}]'); title('Reconstructed (normalised units)');
    % - general stuff
    mySgTitle=sprintf('interpolation: %s - filter: %s',myInterpolation,myFilter);
    if (exist("myTitle","var")), mySgTitle=sprintf("%s - %s",myTitle,mySgTitle); end
    sgtitle(mySgTitle);
end

function ShowSinogramProfiles(profiles,angles)
    ff=figure();
    PlotSinogramProfiles(profiles,angles);
    grid(); colorbar;
    xlabel('Angle [degs]'); ylabel("z [mm]"); title('sinogram');
end

function PlotSinogramProfiles(profiles,myPar)
    PlotSinogram(myPar,profiles(:,1),profiles(:,2:end));
end

function PlotSinogram(xs,ys,data)
    imagesc('XData',xs,'YData',ys,'CData',data);
end
