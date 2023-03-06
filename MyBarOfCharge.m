% {}~

%% include libraries
% - include Matlab libraries
pathToLibrary="externals\MatLabTools";
addpath(genpath(pathToLibrary));
% - include lib libraries
pathToLibrary="lib";
addpath(genpath(pathToLibrary));

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
angle=30;
RM=Rot2D(-angle);
% % - in normalised phase space
% lNorm=true;
% MyPointsStartTransported=AdvanceMyPoints(MyPointsStart,RM,lNorm);
% - in physical phase space
lNorm=false;
MyPointsStartTransported=AdvanceMyPoints(MyPointsStart,RM,lNorm);

%% acquire transfer matrices
% ScanCase="2022-03-13\U1-008A-QUE_C270_secondoGiro_CAM";
% ScanCase="2022-03-13\U1-014A-QUE_C270_secondoGiro_CAM";
% ScanCase="2022-03-13\U1-018A-QUE_C270_secondoGiro_CAM";
% ScanCase="2022-03-13\U2-006A-QUE_C270_secondoGiro_CAM";
% ScanCase="2022-03-13\U2-010A-QUE_C270_secondoGiro_CAM";
ScanCase="2022-03-13\U2-016A-QUE_C270_secondoGiro_CAM";
reMatFileName=strcat("D:\emittanzeHEBT\externals\optics\HEBT\",ScanCase,"_ReMat.tfs");
fprintf('parsing file %s ...\n',reMatFileName);
clear rMatrix; rMatrix = readmatrix(reMatFileName,'HeaderLines',1,'Delimiter',',','FileType','text');
clear TMs; TMs=GetTransMatrix(rMatrix,"HOR","SCAN");

%% acquire profiles
measPath="S:\Area Ricerca\EMITTANZE SUMMARY\EMITTANZE SUMMARY";
iCurr2mon=[-2 -1]; % 1=CAM, 2=DDS; iCAM=iCurr+iCurr2mon(1); iDDS=iCurr+iCurr2mon(2);
iScanSetUps=1;
iMon=1;
% run("D:\emittanzeHEBT\2022-03-13\SetMeUp_U1p008ApQUE_C270_secondoGiro.m");
% run("D:\emittanzeHEBT\2022-03-13\SetMeUp_U1p014ApQUE_C270_secondoGiro.m");
% run("D:\emittanzeHEBT\2022-03-13\SetMeUp_U1p018ApQUE_C270_secondoGiro.m");
% run("D:\emittanzeHEBT\2022-03-13\SetMeUp_U2p006ApQUE_C270_secondoGiro.m");
% run("D:\emittanzeHEBT\2022-03-13\SetMeUp_U2p010ApQUE_C270_secondoGiro.m");
run("D:\emittanzeHEBT\2022-03-13\SetMeUp_U2p016ApQUE_C270_secondoGiro.m");
[measProfiles,CyCodesProf,CyProgsProf]=ParseDDSProfiles(sprintf("%s\\%s\\%s\\profiles\\*_profiles.txt",measPath,path(iScanSetUps),CAMpaths(1)),"CAM");
[BARsMeas,FWxMsMeas,INTsMeas]=StatDistributions(measProfiles);

%% rotate Bar
lNorm=true;
angles=-180:10:180;
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
ShowScans(angles,"\theta [deg]",alphaFit,betaFit,emiGfit,phisFit,sigZ,sigZP,"","Rotating Bar");
ComparePropagationVsFit(angles,"\theta [deg]",alphaFit,betaFit,emiGfit,phisFit,alphaPrp,betaPrp,emiGprp,phisPrp,sigZ,sigZP,"Rotating Bar");

%% try to reproduce scan
angles=20:5:70; % [degs]
nAngles=length(angles);
clear RMs; RMs=Rot2D(-angles);

nPoints=10000;
BBs=6.75*1E-3;    % [m]
HHs=0.025*1E-3;  % [rad]

rr=150; % [m]
xb=linspace(-rr,rr,2*round(rr)+1); % [mm]
yb=linspace(-rr,rr,2*round(rr)+1); % [mrad]

clear alphaStart betaStart emiGstart; alphaStart=NaN(nAngles,length(BBs)); betaStart=NaN(nAngles,length(BBs)); emiGstart=NaN(nAngles,length(BBs));
clear alphaFit betaFit emiGfit phisFit; alphaFit=NaN(size(TMs,3),nAngles,length(BBs)); betaFit=NaN(size(TMs,3),nAngles,length(BBs)); emiGfit=NaN(size(TMs,3),nAngles,length(BBs)); phisFit=NaN(size(TMs,3),nAngles,length(BBs));
clear alphaPrp betaPrp emiGprp phisPrp; alphaPrp=NaN(size(TMs,3),nAngles,length(BBs)); betaPrp=NaN(size(TMs,3),nAngles,length(BBs)); emiGprp=NaN(size(TMs,3),nAngles,length(BBs)); phisPrp=NaN(size(TMs,3),nAngles,length(BBs));
clear sigZ sigZP; sigZ=NaN(size(TMs,3),nAngles,length(BBs)); sigZP=NaN(size(TMs,3),nAngles,length(BBs));
clear profiles; profiles=NaN(length(xb)-1,1+size(TMs,3),nAngles,length(BBs));

% lFirst=true;
for iDim=1:length(BBs)
    MyPointsStart=SampleRect(nPoints,BBs(iDim),HHs(iDim));
    % rotate in physical phase space
    normAngle=-6; % [degs]
    clear RMN; RMN=Rot2D(-normAngle);
    MyPointsStart=AdvanceMyPoints(MyPointsStart,RMN,false,false); % lNorm,ldebug
    % rotate in normalised phase space
    for iRotMat=1:nAngles
        clear RotatedPoints; RotatedPoints=AdvanceMyPoints(MyPointsStart,RMs(:,:,iRotMat),true,false); % lNorm,ldebug
        [alphaStart(iRotMat,iDim),betaStart(iRotMat,iDim),emiGstart(iRotMat,iDim)]=GetOpticsFromSigmaMatrix(RotatedPoints);   % ellypse orientation
        emiGstart(iRotMat,iDim)=max(GetSinglePartEmittance(RotatedPoints,alphaStart(iRotMat,iDim),betaStart(iRotMat,iDim)));  % max emittance

        [alphaFit(:,iRotMat,iDim),betaFit(:,iRotMat,iDim),emiGfit(:,iRotMat,iDim),phisFit(:,iRotMat,iDim),...
         alphaPrp(:,iRotMat,iDim),betaPrp(:,iRotMat,iDim),emiGprp(:,iRotMat,iDim),phisPrp(:,iRotMat,iDim),...
         sigZ(:,iRotMat,iDim),sigZP(:,iRotMat,iDim),profiles(:,:,iRotMat,iDim)]=...
            ScanBar(TMs(1:2,1:2,:),RotatedPoints,xb,yb,betaStart(iRotMat,iDim),alphaStart(iRotMat,iDim),false);%,sprintf("Rotation matrix: %d degs",angles(iRotMat))); % lNorm
    end
end

% Xs=1:size(TMs,3); xLab="ID []"; myLeg=compose("\\theta=%g degs",angles);
% for iDim=1:length(BBs)
%     caseName=sprintf("scan %s - Rotating Bar - bb=%g m; hh=%g rad;",ScanCase,BBs(iDim),HHs(iDim));
%     Show3Ddistributions(profiles(:,:,:,iDim),Xs,xLab,myLeg,caseName);
%     ShowScans(Xs,xLab,alphaFit(:,:,iDim),betaFit(:,:,iDim),emiGfit(:,:,iDim),phisFit(:,:,iDim),sigZ(:,:,iDim),sigZP(:,:,iDim),myLeg,caseName);
% %     ComparePropagationVsFit(Xs,xLab,alphaFit(:,:,iDim),betaFit(:,:,iDim),emiGfit(:,:,iDim),phisFit(:,:,iDim),...
% %                                     alphaPrp(:,:,iDim),betaPrp(:,:,iDim),emiGprp(:,:,iDim),phisPrp(:,:,iDim),...
% %                                     sigZ(:,:,iDim),sigZP(:,:,iDim),strcat(ScanCase," - ",myLeg));
% end

%% find FWHM of simulated profiles
BARsSimul=NaN(size(TMs,3),nAngles,length(BBs)); FWxMsSimul=NaN(size(TMs,3),nAngles,length(BBs)); INTsSimul=NaN(size(TMs,3),nAngles,length(BBs));
for iDim=1:length(BBs)
    for iRotMat=1:nAngles
        [BARsSimul(:,iRotMat,iDim),FWxMsSimul(:,iRotMat,iDim),INTsSimul(:,iRotMat,iDim)]=StatDistributions(profiles(:,:,iRotMat,iDim));
    end
end

%% find FWHM of measured profiles
iPlane=1;
iMin=indices(1+iMon,1,iScanSetUps)-indices(1+iMon,1,iScanSetUps)+1;
iMax=indices(1+iMon,2,iScanSetUps)-indices(1+iMon,1,iScanSetUps)+1;
xLab="ID []"; myLeg=compose("\\theta=%g degs",angles);
for iDim=1:length(BBs)
    figure(); cm=colormap(parula(nAngles));
    for ii=1:nAngles
        if ( angles(ii)<0 )
            plot(iMin:iMax,FWxMsSimul(:,ii,iDim),"^-","Color",cm(ii,:));
        elseif ( angles(ii)==0 )
            plot(iMin:iMax,FWxMsSimul(:,ii,iDim),"o-","Color",cm(ii,:),"MarkerSize",10);
        else
            plot(iMin:iMax,FWxMsSimul(:,ii,iDim),"v-","Color",cm(ii,:));
        end
        hold on;
    end
    plot(iMin:iMax,FWxMsMeas(indices(1+iMon,1,iScanSetUps):indices(1+iMon,2,iScanSetUps),iPlane),"k*-");
    grid(); xlabel(xLab); ylabel("FWHM [mm]"); legend([myLeg "measurement"],"Location","best");
    title(sprintf("%s - bb=%g; hh=%g",ScanCase,BBs(iDim),HHs(iDim)));
    FigFileOut=sprintf("%s_BB%g_HH%g_tMin%g_tMax%g_phys-6_norm.fig",ScanCase,BBs(iDim),HHs(iDim),min(angles),max(angles));
    savefig(FigFileOut);
    fprintf("...saving to file %s ...\n",FigFileOut);
end

%% compare simulation and measurements
iPlane=1; iFitRange=1;
myProfIDs=indices(1+iMon,1,iScanSetUps):indices(1+iMon,2,iScanSetUps);
allCurrIDs=indices(1,1,iScanSetUps):indices(1,2,iScanSetUps);
% iMin=1; iMax=fitIndices(iMon,2,iPlane,iFitRange,iScanSetUps)-fitIndices(iMon,1,iPlane,iFitRange,iScanSetUps)+1;
nProfiles=length(myProfIDs);
[nRows,nCols]=GetNrowsNcols(nProfiles+1);
ff=figure();
ff.Position(1:2)=[0 0]; % figure at lower-left corner of screen
ff.Position(3)=ff.Position(3)*nCols/2; % larger figure
ff.Position(4)=ff.Position(4)*nRows/2; % larger figure
cm=colormap(parula(nAngles));
t = tiledlayout(nRows,nCols,'TileSpacing','Compact','Padding','Compact'); % minimise whitespace around plots
iShift=0;
for iProf=1:nProfiles
    nexttile;
    % plot data
    % myProfID=fitIndices(iMon,1,iPlane,iFitRange,iScanSetUps)+(iProf-1);
    myCurrID=find(allCurrIDs==myProfIDs(iProf)+(indices(1,1,iScanSetUps)-indices(1+iMon,1,iScanSetUps)))+iShift;
    if ( 1<=myCurrID && myCurrID<=size(profiles,2)-1)
        for ii=1:nAngles
            [myMax,myMaxID]=max(profiles(:,1+myCurrID,ii));
            plot(profiles(:,1,ii),profiles(:,1+myCurrID,ii)/myMax,".-","Color",cm(ii,:)); hold on;
        end
    end
%     [myMax,myMaxID]=max(measProfiles(:,1+myProfIDs(iProf),iPlane));
%     plot(measProfiles(:,1,iPlane)-measProfiles(myMaxID,1,iPlane),measProfiles(:,1+myProfIDs(iProf),iPlane)/myMax,"k.-");
    [myMax,myMaxID]=max(measProfiles(:,1+myProfIDs(iProf),iPlane));
    plot(measProfiles(:,1,iPlane)-BARsMeas(myProfIDs(iProf),iPlane),measProfiles(:,1+myProfIDs(iProf),iPlane)/myMax,"k.-");
    title(sprintf("profile ID: CAM=%d; ID=%g",myProfIDs(iProf),iProf));
    xlabel("position [mm]");
    ylabel("profile normalised to max");
    grid();
    if (~isnan(FWxMsMeas(myProfIDs(iProf),iPlane)))
        xlim([-FWxMsMeas(myProfIDs(iProf),iPlane) FWxMsMeas(myProfIDs(iProf),iPlane)]);
    end
end
% legend
nexttile;
for ii=1:nAngles
    plot(NaN(),NaN(),".-","Color",cm(ii,:)); hold on;
end
plot(NaN(),NaN(),"k.-");
legend([myLeg "measured"],"Location","best");
sgtitle(ScanCase);

%% functions

function [alpha,beta,emiG]=ReconstructOptics(myCoords,algo)
    if (~exist("algo","var")), algo="inscribed"; end
    switch upper(algo)
        case "INSCRIBED"
            alpha=0.0;
            xMax=max(myCoords(:,1));
            xpMax=max(myCoords(:,2));
            beta=xMax/xpMax;
            emiG=xMax*xpMax;
        case "CIRCUMSCRIBED"
            [alpha,beta,emiG]=ReconstructOptics(myCoords,"INSCRIBED");
            emiG=2*emiG;
        otherwise
            error("...reconstruction algorithm %s NOT recognised.",algo);
    end
end

function [zOut]=RotateMe(zIn,angle)
    TT=[cosd(angle) -sind(angle); sind(angle) cosd(angle)];
    zOut=TT*zIn';
    zOut=zOut';
end

function ShowMe(MyPoints,MyContours)
    figure();
    
    % points
    plot(MyPoints(:,1),MyPoints(:,2),"k.");
    
    % contours
    hold on;
    for ii=1:size(MyContours,3)
        if (ii>1), hold on; end
        plot(MyContours(:,1,ii),MyContours(:,2,ii),".-");
    end
    
    grid(); % axis equal;
end
