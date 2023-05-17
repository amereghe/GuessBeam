% {}~

%% include libraries
% % - include Matlab libraries
% pathToLibrary="../MatLabTools";
% addpath(genpath(pathToLibrary));
% % - include lib libraries
% pathToLibrary="lib";
% addpath(genpath(pathToLibrary));

%% 
clear all; % clean RAM
close all; % close all figures

%% user settings
whatScan="MUX";
planeShow="X";

% - starting optics functions: V.Lante, "XPR optics design"
BETXstartNom=8.184; ALFXstartNom=-0.443;
BETYstartNom=4.526; ALFYstartNom=-1.985;
emiGeoNom=1E-6; % [m rad]

% - path to files
% iFileNameOpticsSummary="../optics/HEBT/output_p030_MUX_X3-011B-VWN_free_10p0degs/x3_011b_vwn_summary_optics.tfs";
% iFileNameRMatrixSummary="../optics/HEBT/output_p030_MUX_X3-011B-VWN_free_10p0degs/x3_011b_vwn_summary_rmatrix.tfs";
% iFileNameOpticsSummary="../optics/HEBT/output_p030_MUY_X3-011B-VWN_free_10p0degs/x3_011b_vwn_summary_optics.tfs";
% iFileNameRMatrixSummary="../optics/HEBT/output_p030_MUY_X3-011B-VWN_free_10p0degs/x3_011b_vwn_summary_rmatrix.tfs";
% iFileNameOpticsSummary="../optics/HEBT/output_p030_MUX_X3-011B-VWN_constrained_10p0degs/x3_011b_vwn_summary_optics.tfs";
% iFileNameRMatrixSummary="../optics/HEBT/output_p030_MUX_X3-011B-VWN_constrained_10p0degs/x3_011b_vwn_summary_rmatrix.tfs";
% iFileNameOpticsSummary="../optics/HEBT/output_p030_MUX_X3-011B-VWN_constrained_02p5degs/x3_011b_vwn_summary_optics.tfs";
% iFileNameRMatrixSummary="../optics/HEBT/output_p030_MUX_X3-011B-VWN_constrained_02p5degs/x3_011b_vwn_summary_rmatrix.tfs";
iFileNameOpticsSummary="../optics/HEBT/output_p320_MUX_X3-011B-VWN_free_10p0degs/x3_011b_vwn_summary_optics.tfs";
iFileNameRMatrixSummary="../optics/HEBT/output_p320_MUX_X3-011B-VWN_free_10p0degs/x3_011b_vwn_summary_rmatrix.tfs";

% - 2D histograms (Poincaré phase plots)
xb=linspace(-64,64,128+1); % [mm]
yb=linspace(-3,3,100+1); % [mrad]

% - beam sampling
nParticles=10000;

% distribution type per plane
iDist="bar";
% is beam distribution sampled in normal coordinates or not?
iNorm=true;
% 
betaStartUsed=BETXstartNom;
alfaStartUsed=ALFXstartNom;
gammaStartUsed=(1+alfaStartUsed.^2)./betaStartUsed;
emiGeoUsed=emiGeoNom;

% - bar of charge dimensions
bb=1;     % [sqrt(m)] or [m]
hh=0.01;  % [sqrt(m)] or [rad]
barIsNormal=false; % bar of charge sampled in normal coordinates

% - Gaussian sampling
GaussIsNormal=true;

%% do the job

% - crunch optics file
scannedOptics=ParseTfsTable(iFileNameOpticsSummary,'optics',"SCAN");
[ colNames, colUnits, colFacts, mapping, readFormat ] = ...
        GetColumnsAndMappingTFS("optics","SCAN");
iID=mapping(strcmp(colNames,'ID'));
iMU=mapping(strcmp(colNames,whatScan));
iMUX=mapping(strcmp(colNames,"MUX"));
iMUY=mapping(strcmp(colNames,"MUY"));
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
MUs=rad2deg(scannedOptics{iMU}*2*pi);
[uMUs,iuMUs,~]=unique(MUs(sortedIDs));
nActualPoints=length(MUs(sortedIDs(iuMUs)));

% - prepare variables
mu0=rad2deg(scannedOptics{iMU}(zeroPos(1))*2*pi);
betaObsNom=NaN(length(planeShow),1); alfaObsNom=NaN(length(planeShow),1);
betaStartNom=NaN(length(planeShow),1); alfaStartNom=NaN(length(planeShow),1);
emigNom=NaN(length(planeShow),1);
for iPlane=1:length(planeShow)
    % - get nominal optics at start of lattice
    if (strcmpi(planeShow(iPlane),"X"))
        betaStartNom(iPlane)=BETXstartNom; alfaStartNom(iPlane)=ALFXstartNom;
    elseif (strcmpi(planeShow(iPlane),"Y"))
        betaStartNom(iPlane)=BETYstartNom; alfaStartNom(iPlane)=ALFYstartNom;
    end
    emigNom(iPlane)=emiGeoNom; % [m rad]
    % - get nominal optics at measurement point
    betCur=sprintf("BET%s",planeShow(iPlane)); alfCur=sprintf("ALF%s",planeShow(iPlane));
    betaObsNom(iPlane)=scannedOptics{mapping(strcmp(colNames,betCur))}(zeroPos(1));
    alfaObsNom(iPlane)=scannedOptics{mapping(strcmp(colNames,alfCur))}(zeroPos(1));
end
gammaStartNom=(1+alfaStartNom.^2)./betaStartNom;
gammaObsNom=(1+alfaObsNom.^2)./betaObsNom;

% - crunch response matrix file
scannedRMatrices=ParseTfsTable(iFileNameRMatrixSummary,'RMatrix',"SCAN");
[ colNames, colUnits, colFacts, mapping, readFormat ] = ...
        GetColumnsAndMappingTFS("RMatrix","SCAN");
% - get RMs
readRMs=missing();
for iPlane=1:length(planeShow)
    if (strcmpi(planeShow(iPlane),"X"))
        scanPlane="HOR";
    elseif (strcmpi(planeShow(iPlane),"Y"))
        scanPlane="VER";
    end
    readRMs=ExpandMat(readRMs,GetTransMatrix(cell2mat(scannedRMatrices),scanPlane,"SCAN"));
end

%% define beam

MyPointsStart=missing();
MyContoursStartDist=missing(); MyContoursStartFit=missing(); MyContoursStartTheory=missing();

for iPlane=1:length(planeShow)
    
    switch upper(iDist(iPlane))
        case "BAR"
            % - sample
            MyPointsStart=ExpandMat(MyPointsStart,BeamSample_Rect(nParticles,bb(iPlane),hh(iPlane)));
            % - contour
            MyContoursStartDist=ExpandMat(MyContoursStartDist,GenPointsAlongRectangle(bb(iPlane),hh(iPlane)));
            if (iNorm(iPlane))
                MyPointsStart(:,:,iPlane)=Norm2Phys(MyPointsStart(:,:,iPlane),betaStartUsed(iPlane),alfaStartUsed(iPlane),emiGeoUsed(iPlane));
                MyContoursStartDist(:,:,iPlane)=Norm2Phys(MyContoursStartDist(:,:,iPlane),betaStartUsed(iPlane),alfaStartUsed(iPlane),emiGeoUsed(iPlane));
            end
        case "GAUSS"
            % - sample
            MyPointsStart=ExpandMat(MyPointsStart,BeamSample_Gauss(nParticles,alfaStartUsed(iPlane),betaStartUsed(iPlane),emiGeoUsed(iPlane)));
            % - contour
            MyContoursStartDist=ExpandMat(MyContoursStartDist,GenPointsAlongEllypse(alfaStartUsed(iPlane),betaStartUsed(iPlane),emiGeoUsed(iPlane)));
    end

    % contouring ellypse
    [alphaStart,betaStart,emiGstart]=GetOpticsFromSigmaMatrix(MyPointsStart(:,:,iPlane));   % ellypse orientation
    emiGstart=max(GetSinglePartEmittance(MyPointsStart(:,:,iPlane),alphaStart,betaStart));  % max 
    MyContoursStartFit=ExpandMat(MyContoursStartFit,GenPointsAlongEllypse(alphaStart,betaStart,emiGstart));
    
    % theoretical ellypse
    MyContoursStartTheory=ExpandMat(MyContoursStartTheory,GenPointsAlongEllypse(alfaStartNom(iPlane),betaStartNom(iPlane),emigNom(iPlane)));

    % show beam sampled
    MyContours=missing();
    MyContours=ExpandMat(MyContours,MyContoursStartDist(:,:,iPlane));
    MyContours=ExpandMat(MyContours,MyContoursStartFit(:,:,iPlane));
    MyContours=ExpandMat(MyContours,MyContoursStartTheory(:,:,iPlane));
    MyContoursLabels=["Beam Generation" "Statistics on Generated Beam" "Nominal Optics"];
    ShowCoordsPhysNormScatter(MyPointsStart(:,:,iPlane),betaStartUsed(iPlane),alfaStartUsed(iPlane),emiGeoUsed(iPlane),sprintf("starting beam - plane %s",planeShow(iPlane)),true,MyContours,MyContoursLabels);
    
    % show 2D histogram
    ShowCoordsPhysNorm(MyPointsStart(:,:,iPlane),betaStartUsed(iPlane),alfaStartUsed(iPlane),emiGeoUsed(iPlane),sprintf("starting beam - plane %s",planeShow(iPlane)));
end

%% scan beam,
lNorm=false; % lNorm=true for normalised units
l3D=true;
angles=MUs(sortedIDs(iuMUs));
RMs=missing(); profiles=missing();
for iPlane=1:length(planeShow)
    % - response matrices
    RMs=ExpandMat(RMs,readRMs(1:2,1:2,sortedIDs(iuMUs),iPlane));
    % - track beam
    [alphaFit,betaFit,emiGfit,phisFit,alphaPrp,betaPrp,emiGprp,phisPrp,sigZ,sigZP,tmpProfiles]=...
        ScanBar(RMs(:,:,:,iPlane),MyPointsStart(:,:,iPlane),xb,yb,betaStartUsed(iPlane),alfaStartUsed(iPlane),lNorm,sprintf("Tracked beam - plane %s",planeShow(iPlane)));
    profiles=ExpandMat(profiles,tmpProfiles);
    Show3Ddistributions(profiles(:,:,iPlane),angles,"\theta [deg]","",sprintf("Tracked beam - plane %s",planeShow(iPlane)));
    % ShowScans(angles,"\theta [deg]",alphaFit,betaFit,emiGfit,phisFit,sigZ,sigZP,"",sprintf("Tracked beam - plane %s",planeShow(iPlane)));
    % ComparePropagationVsFit(angles,"\theta [deg]",alphaFit,betaFit,emiGfit,phisFit,alphaPrp,betaPrp,emiGprp,phisPrp,sigZ,sigZP,sprintf("Tracked beam - plane %s",planeShow(iPlane)));

    % show simulated profiles as a sinogram
    ShowSinogramProfiles(profiles(:,:,iPlane),angles,l3D,sprintf("Tracked beam - plane %s",planeShow(iPlane)));
end

%% transport beam through nominal optics
RM=RMs(:,:,mu0==angles,:);
lNorm=false; % lNorm=true for normalised units
MyPointsTransportedNominal=missing();
for iPlane=1:length(planeShow)
    MyPointsTransportedNominal=ExpandMat(MyPointsTransportedNominal,AdvanceMyPoints(MyPointsStart(:,:,iPlane),RM(:,:,:,iPlane),lNorm));
    % show 2D histogram
    ShowCoordsPhysNorm(MyPointsTransportedNominal(:,:,iPlane),betaObsNom(iPlane),alfaObsNom(iPlane),emigNom(iPlane),sprintf("Tracked beam through nominal optics - plane %s",planeShow(iPlane)));
end

%% try to get inverse radon transform
% - select range of points
myRange=1:nActualPoints;
% - manipulate angles:
%   mu=-Radon_angle
%   subtract mu0 if you want to compare to transported beam
myAngles=-(angles(myRange));
% - manipulate profiles:
%   selected range
myProfiles=profiles(:,[1 myRange+1],:);
%   normalised coordinates
for iPlane=1:length(planeShow)
    myProfiles(:,1,iPlane)=myProfiles(:,1,iPlane)/(sqrt(betaObsNom(iPlane)*emigNom(iPlane))*1E3);
    ShowSinogramIR(myProfiles(:,:,iPlane),myAngles);
end

%% stop
return

%% scan angular interval
Delta=180;
for iPlane=1:length(planeShow)
    iMin=0;
    for iMax=find(diff(angles>angles(1)+Delta)):length(angles)
        iMin=iMin+1;
        myRange=iMin:iMax;
        myTitle=sprintf("\\theta[degs]=[%g:%g] - \\Delta\\theta[degs]=%g - plane %s",angles(iMin),angles(iMax),angles(iMax)-angles(iMin),planeShow(iPlane));
        ShowSinogramIR(myProfiles(:,[1 myRange+1],iPlane),myAngles(myRange),myTitle);
        pause(1);
    end
end

%% local functions
function ShowCoordsPhysNormScatter(MyPoints,beta,alfa,emiGeo,myTitle,lPhys,MyContours,MyContoursLabel)
    if (~exist("lPhys","var")), lPhys=true; end
    MyContoursPhys=missing(); MyContoursNorm=missing();
    if (lPhys)
        MyPointsPhys=MyPoints;
        MyPointsNorm=Phys2Norm(MyPointsPhys,beta,alfa,emiGeo);
        if (exist("MyContours","var")), MyContoursPhys=MyContours; end
    else
        MyPointsNorm=MyPoints;
        MyPointsPhys=Norm2Phys(MyPointsNorm,beta,alfa,emiGeo);
        if (exist("MyContours","var")), MyContoursNorm=MyContours; end
    end
    
    %% compute histograms
    clear nCountsPhys nCountsXphys nCountsYphys xbPhys ybPhys;
    [nCountsPhys,nCountsXphys,nCountsYphys,xbPhys,ybPhys]=Get2dHistograms(MyPointsPhys(:,1),MyPointsPhys(:,2));
    clear nCountsNorm nCountsXnorm nCountsYnorm xbNorm ybNorm;
    [nCountsNorm,nCountsXnorm,nCountsYnorm,xbNorm,ybNorm]=Get2dHistograms(MyPointsNorm(:,1),MyPointsNorm(:,2));
    
    %% actually show
    figure();
    [axM,axX,axY]=Plot2DHistograms(MyPointsPhys,nCountsXphys,nCountsYphys,xbPhys,ybPhys,"z [m]","zp [rad]",MyContoursPhys,false);
    sgtitle(sprintf("%s - Physical units",myTitle));
    if (~all(ismissing(MyContoursPhys),"all") & exist("MyContoursLabel","var")), legend(axM,["Beam" MyContoursLabel],"Location","best"); end
    figure();
    [axM,axX,axY]=Plot2DHistograms(MyPointsNorm,nCountsXnorm,nCountsYnorm,xbNorm,ybNorm,"z []","zp []",MyContoursNorm,false);
    sgtitle(sprintf("%s - Normalised units",myTitle));
    if (~all(ismissing(MyContoursNorm),"all") & exist("MyContoursLabel","var")), legend(axM,["Beam" MyContoursLabel],"Location","best"); end
end

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
    PlotSinogramProfiles(profiles,angles,l3D);
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

function ShowSinogramProfiles(profiles,angles,l3D,myTitle)
    if (~exist("l3D","var")), l3D=true; end
    ff=figure();
    PlotSinogramProfiles(profiles,angles,l3D);
    grid(); colorbar;
    xlabel('Angle [degs]'); ylabel("z [mm]");
    if (exist("myTitle","var"))
        title(myTitle);
    else
        title('sinogram');
    end
end

function PlotSinogramProfiles(profiles,myPar,l3D)
    PlotSinogram(myPar,profiles(:,1),profiles(:,2:end),l3D);
end

function PlotSinogram(xs,ys,data,l3D)
    if (~exist("l3D","var")), l3D=true; end
    if (l3D)
        % 3D plot
        surf(xs,ys,data,"EdgeColor","none");
        view(2); % starts appearing as a 2D plot
    else
        % 2D plot
        imagesc('XData',xs,'YData',ys,'CData',data);
        colorbar;
    end
    grid();
end
