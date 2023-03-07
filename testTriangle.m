%% generate distribution
aa=3;
bb=5;
cc=8;
dd=13;
pp=[aa bb cc dd];
myGrid=101;
shape="trapezoid";

%% sample according to distribution
nPoints=100000;

% - actual pdf
xx=linspace(min(pp),max(pp),myGrid);
yy=EvalMe(xx,pp,shape);

% - cumulative pdf
YCUM=cumtrapz(yy); YCUM=YCUM/YCUM(end);
%   debug cumulative pdf
% figure();
% plot(xx,yy,".-");
% yyaxis right;
% plot(xx,YCUM,".-");
% yyaxis left;
% figure();
% plot(YCUM,xx,".-");

% - actually sample
rng('default'); rng(1);
xSample=interp1(YCUM,xx,rand(nPoints,1),"spline");

% - get histogram
[myHist,xb]=histcounts(xSample,100,"Normalization","pdf");
edges=xb(2:end)-0.5*(xb(2)-xb(1)); 

%%  show everything
figure();
% - histogram
bar(edges,myHist,1);
% - original distribution
hold on; plot(xx,yy,"r.-");
% - analytical stats
yl=ylim(); yE=mean(yl);
lGeneral=true; muG=AnalyticalMean(pp,shape,lGeneral); sig2G=AnalyticalSigma2(pp,shape,lGeneral);
hold on; errorbar(muG,yE,sqrt(sig2G),"horizontal","ms");
lGeneral=false; muP=AnalyticalMean(pp,shape,lGeneral); sig2P=AnalyticalSigma2(pp,shape,lGeneral);
hold on; errorbar(muP,yE*0.9,sqrt(sig2P),"horizontal","c*");
muS=mean(xSample); sigS=std(xSample);
hold on; errorbar(muS,yE*0.8,sigS,"horizontal","g*");
% - general stuff
legend("histogram sampled population","analytical profile",...
    "analytical stat (general)","analytical stat (spec)","stat on population","Location","best");
grid on;

%% local functions

function hh=MyHeight(pp,shape)
    if (~exist("shape","var")), shape="rtriangle"; end
    switch upper(shape)
        case "TRIANGLE"
            % A=(aa,0); B=(bb,hh); C=(cc,0);
            % pp(1)=aa; pp(2)=bb; pp(3)=cc;
            hh=2/(pp(3)-pp(1));
        case "TRAPEZOID"
            % A=(aa,0); B=(bb,hh); C=(cc,hh); D=(dd,0);
            % pp(1)=aa; pp(2)=bb; pp(3)=cc; pp(4)=dd;
            hh=2/(pp(3)-pp(2)+pp(4)-pp(1));
        case {"RTRIANGLE","LTRIANGLE"}
            % R: A=(aa,0); B=(bb,hh);
            % L: A=(aa,hh); B=(bb,0);
            % pp(1)=aa; pp(2)=bb;
            hh=2/(pp(2)-pp(1));
        case "RECTANGLE"
            % pp(1)=aa; pp(2)=bb;
            hh=1/(pp(2)-pp(1));
        otherwise
            error("...unknown shape %s.",shape);
    end
end

function YY=EvalMe(XX,pp,shape)
    if (~exist("shape","var")), shape="rtriangle"; end
    hh=MyHeight(pp,shape);
    switch upper(shape)
        case "TRIANGLE"
            % A=(aa,0); B=(bb,hh); C=(cc,0);
            % pp(1)=aa; pp(2)=bb; pp(3)=cc;
            YY=EvalMeBasicShapes(XX,pp(1:2),hh,"RTRIANGLE")+...
               EvalMeBasicShapes(XX,pp(2:3),hh,"LTRIANGLE");
            indices=(XX==pp(2));
            if ( any(indices) )
                YY(indices)=YY(indices)/2;
            end
        case "TRAPEZOID"
            % A=(aa,0); B=(bb,hh); C=(cc,hh); D=(dd,0);
            % pp(1)=aa; pp(2)=bb; pp(3)=cc; pp(4)=dd;
            YY=EvalMeBasicShapes(XX,pp(1:2),hh,"RTRIANGLE")+...
               EvalMeBasicShapes(XX,pp(2:3),hh,"RECTANGLE")+...
               EvalMeBasicShapes(XX,pp(3:4),hh,"LTRIANGLE");
            for ip=2:3
                indices=(XX==pp(ip));
                if ( any(indices) )
                    YY(indices)=YY(indices)/2;
                end
            end
        otherwise
            YY=EvalMeBasicShapes(XX,pp,hh,shape);
    end
end
function YY=EvalMeBasicShapes(XX,pp,hh,shape)
    if (~exist("shape","var")), shape="rtriangle"; end
    YY=zeros(size(XX));
    switch upper(shape)
        case "RTRIANGLE"
            % A=(aa,0); B=(bb,hh);
            % pp(1)=aa; pp(2)=bb;
            indices=(pp(1)<=XX & XX<=pp(2));
            YY(indices)=hh*(XX(indices)-pp(1))/(pp(2)-pp(1));
        case "LTRIANGLE"
            % A=(aa,hh); B=(bb,0);
            % pp(1)=aa; pp(2)=bb;
            indices=(pp(1)<=XX & XX<=pp(2));
            YY(indices)=hh*(XX(indices)-pp(2))/(pp(1)-pp(2));
        case "RECTANGLE"
            % A=(aa,hh); B=(bb,hh);
            % pp(1)=aa; pp(2)=bb;
            indices=(pp(1)<=XX & XX<=pp(2));
            YY(indices)=hh;
        otherwise
            error("...unknown shape %s.",shape);
    end
end

function mu=AnalyticalMean(pp,shape,lGeneral)
    if (~exist("shape","var")), shape="rtriangle"; end
    if (~exist("lGeneral","var")), lGeneral=true; end
    hh=MyHeight(pp,shape);
    switch upper(shape)
        case "TRIANGLE"
            % A=(aa,0); B=(bb,hh); C=(cc,0);
            % pp(1)=aa; pp(2)=bb; pp(3)=cc;
            mu=AnalyticalMeanBasicShapes(pp(1:2),"RTRIANGLE",hh,lGeneral)+...
               AnalyticalMeanBasicShapes(pp(2:3),"LTRIANGLE",hh,lGeneral);
        case "TRAPEZOID"
            % A=(aa,0); B=(bb,hh); C=(cc,hh); D=(dd,0);
            % pp(1)=aa; pp(2)=bb; pp(3)=cc; pp(4)=dd;
            mu=AnalyticalMeanBasicShapes(pp(1:2),"RTRIANGLE",hh,lGeneral)+...
               AnalyticalMeanBasicShapes(pp(2:3),"RECTANGLE",hh,lGeneral)+...
               AnalyticalMeanBasicShapes(pp(3:4),"LTRIANGLE",hh,lGeneral);
        otherwise
            mu=AnalyticalMeanBasicShapes(pp,shape,hh,lGeneral);
    end
end
function mu=AnalyticalMeanBasicShapes(pp,shape,hh,lGeneral)
    if (~exist("shape","var")), shape="rtriangle"; end
    if (~exist("lGeneral","var")), lGeneral=true; end
    switch upper(shape)
        case "RTRIANGLE"
            % A=(aa,0); B=(bb,hh);
            % pp(1)=aa; pp(2)=bb;
            if (lGeneral) 
                mu=hh/6*(pp(2)-pp(1))*(2*pp(2)+pp(1));
            else
                mu=(2*pp(2)+pp(1))/3;
            end
        case "LTRIANGLE"
            % A=(aa,hh); B=(bb,0);
            % pp(1)=aa; pp(2)=bb;
            if (lGeneral) 
                mu=hh/6*(pp(2)-pp(1))*(2*pp(1)+pp(2));
            else
                mu=(2*pp(1)+pp(2))/3;
            end
        case "RECTANGLE"
            % A=(aa,hh); B=(bb,hh);
            % pp(1)=aa; pp(2)=bb;
            if (lGeneral) 
                mu=hh/2*(pp(2)-pp(1))*(pp(2)+pp(1));
            else
                mu=(pp(1)+pp(2))/2;
            end
        otherwise
            error("...unknown shape %s.",shape);
    end
end

function sig2=AnalyticalSigma2(pp,shape,lGeneral)
    if (~exist("shape","var")), shape="rtriangle"; end
    if (~exist("lGeneral","var")), lGeneral=true; end
    hh=MyHeight(pp,shape); mu=AnalyticalMean(pp,shape,lGeneral);
    switch upper(shape)
        case "TRIANGLE"
            % A=(aa,0); B=(bb,hh); C=(cc,0);
            % pp(1)=aa; pp(2)=bb; pp(3)=cc;
            sig2=AnalyticalSigma2BasicShapes(pp(1:2),"RTRIANGLE",hh,mu,lGeneral)+...
                 AnalyticalSigma2BasicShapes(pp(2:3),"LTRIANGLE",hh,mu,lGeneral);
        case "TRAPEZOID"
            % A=(aa,0); B=(bb,hh); C=(cc,hh); D=(dd,0);
            % pp(1)=aa; pp(2)=bb; pp(3)=cc; pp(4)=dd;
            sig2=AnalyticalSigma2BasicShapes(pp(1:2),"RTRIANGLE",hh,mu,lGeneral)+...
                 AnalyticalSigma2BasicShapes(pp(2:3),"RECTANGLE",hh,mu,lGeneral)+...
                 AnalyticalSigma2BasicShapes(pp(3:4),"LTRIANGLE",hh,mu,lGeneral);
        otherwise
            sig2=AnalyticalSigma2BasicShapes(pp,shape,hh,mu,lGeneral);
    end
end
function sig2=AnalyticalSigma2BasicShapes(pp,shape,hh,mu,lGeneral)
    if (~exist("shape","var")), shape="rtriangle"; end
    if (~exist("lGeneral","var")), lGeneral=true; end
    switch upper(shape)
        case "RTRIANGLE"
            % A=(aa,0); B=(bb,hh);
            % pp(1)=aa; pp(2)=bb;
            if (lGeneral) 
                sig2=hh*(pp(2)-pp(1))/12*((pp(1)+pp(2)-2*mu)^2+2*(pp(2)-mu)^2);
            else
                sig2=(pp(2)-pp(1))^2/18;
            end
        case "LTRIANGLE"
            % A=(aa,hh); B=(bb,0);
            % pp(1)=aa; pp(2)=bb;
            if (lGeneral) 
                sig2=hh*(pp(2)-pp(1))/12*((pp(1)+pp(2)-2*mu)^2+2*(pp(1)-mu)^2);
            else
                sig2=(pp(2)-pp(1))^2/18;
            end
        case "RECTANGLE"
            % A=(aa,hh); B=(bb,hh);
            % pp(1)=aa; pp(2)=bb;
            if (lGeneral) 
                sig2=hh*(pp(2)-pp(1))/3*((pp(1)+pp(2)-2*mu)^2-(pp(1)-mu)*(pp(2)-mu));
            else
                sig2=(pp(2)-pp(1))^2/12;
            end
        otherwise
            error("...unknown shape %s.",shape);
    end
end