% {}~

%% include libraries
% - include Matlab libraries
pathToLibrary="externals\MatLabTools";
addpath(genpath(pathToLibrary));
% - include Matlab libraries
pathToLibrary="externals\ExternalMatLabTools";
addpath(genpath(pathToLibrary));
% - include local library with functions
pathToLibrary="lib";
addpath(genpath(pathToLibrary));

%% define bar of charge and ellypse - norm phase space

% bar of charge
bb=1; nX=31;
hh=0.1; nY=11;
xDom=linspace(-bb,bb,nX);
yDom=linspace(-hh,hh,nY); yDom=yDom(2:end-1); % avoid repetitions
% ellypse
nEl=360;

% build bar of charge
MyContour=NaN(2*(length(xDom)+length(yDom))+1,2);
MyContour(1:end-1,1)=[ xDom                     bb*ones(1,length(yDom)) xDom(end:-1:1)          -bb*ones(1,length(yDom)) ]; 
MyContour(1:end-1,2)=[ -hh*ones(1,length(xDom)) yDom                    hh*ones(1,length(xDom)) yDom(end:-1:1)           ];
MyContour(end,:)=MyContour(1,:);

% build external ellypse
MyEllypse=NaN(nEl+1,2);
rr=max(sqrt(MyContour(:,1).^2+MyContour(:,2).^2));
MyEllypse(1:end-1,1)=cosd(linspace(0,360,nEl))*rr;
MyEllypse(1:end-1,2)=sind(linspace(0,360,nEl))*rr;
MyEllypse(end,:)=MyEllypse(1,:);

ShowMe(MyContour,MyEllypse); axis equal;

%% sample points
% sample points
nPoints=10000;
nBinsX=50;
nBinsY=nBinsX;

MyPoints=NaN(nPoints,2);
MyPoints(:,1)=2*bb*rand(nPoints,1)-bb;
MyPoints(:,2)=2*hh*rand(nPoints,1)-hh;

ShowMe(MyContour,MyEllypse,MyPoints); axis equal;

% show histograms
% xb=linspace(-bb,bb,100);
% yb=linspace(-hh,hh,100);
% [nCounts,nCountsX,nCountsY]=Get2dHistograms(MyPoints(:,1),MyPoints(:,2),xb,yb);
[nCounts,nCountsX,nCountsY,xb,yb]=Get2dHistograms(MyPoints(:,1),MyPoints(:,2));
Show2DHistograms(nCounts,nCountsX,nCountsY,xb,yb,"z []","zp []","Normalised phase space");

%% translate to physical space
beta=1; % [m]
alpha=-2; % []
emiG=1E-6; % [pi m rad]
MyContourReal=Norm2Phys(MyContour,beta,alpha,emiG);
MyEllypseReal=Norm2Phys(MyEllypse,beta,alpha,emiG);

% show stuff
ShowMe(MyContourReal,MyEllypseReal);

%% rotate bar of charge
angle=30; % [degs]
MyContourRotated=RotateMe(MyContour,-angle);
ShowMe(MyContourRotated,MyEllypse); axis equal;

%% translate to physical space
MyContourRotatedReal=Norm2Phys(MyContourRotated,beta,alpha,emiG);

% show stuff
ShowMe(MyContourRotatedReal,MyEllypseReal);

%% functions

function [zOut]=RotateMe(zIn,angle)
    TT=[cosd(angle) -sind(angle); sind(angle) cosd(angle)];
    zOut=TT*zIn';
    zOut=zOut';
end

function [zOut]=Phys2Norm(zIn,beta,alpha,emiG)
    TT=[1 0; alpha beta]/sqrt(beta*emiG);
    zOut=TT*zIn';
    zOut=zOut';
end

function [zOut]=Norm2Phys(zIn,beta,alpha,emiG)
    TT=[1 0; -alpha/beta 1/beta]*sqrt(beta*emiG);
    zOut=TT*zIn';
    zOut=zOut';
end

function ShowMe(MyContour,MyEllypse,MyPoints)
    figure();
    plot(MyContour(:,1),MyContour(:,2),"k.-",MyEllypse(:,1),MyEllypse(:,2),"r-");
    grid(); % axis equal;
    if ( exist("MyPoints","var") )
        hold on;
        plot(MyPoints(:,1),MyPoints(:,2),".");
    end
end