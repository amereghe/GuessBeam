function [alphaFit,betaFit,emiGfit,phisFit,alphaPrp,betaPrp,emiGprp,phisPrp,sigZ,sigZP,profiles]=...
    ScanBar(TMs,StartPoints,xb,yb,betaStart,alphaStart,lNorm,myTitle)
    if (~exist("lNorm","var")), lNorm=false; end
    if (~exist("myTitle","var")), myTitle=missing(); end
    nIter=size(TMs,3);
    alphaFit=NaN(nIter,1); betaFit=NaN(nIter,1); emiGfit=NaN(nIter,1); phisFit=NaN(nIter,1);
    alphaPrp=NaN(nIter,1); betaPrp=NaN(nIter,1); emiGprp=NaN(nIter,1); phisPrp=NaN(nIter,1);
    sigZ=NaN(nIter,1); sigZP=NaN(nIter,1);
    
    profiles=NaN(length(xb)-1,1+nIter);
    profiles(:,1)=xb(1:end-1)+(xb(2)-xb(1))/2;

    if (~ismissing(myTitle))
        ff=figure();
        ff.Position(1:2)=[1000 100]; % figure at lower-left corner of screen
        ff.Position(3:4)=ff.Position(3:4)*1.5; % larger figure
    end
    
    for iIter=1:nIter
        clear MyContours; MyContours=missing();
        TransportedPoints=AdvanceMyPoints(StartPoints,TMs(:,:,iIter),lNorm,false);

        % contouring ellypse
        [alphaFit(iIter),betaFit(iIter),emiGfit(iIter),sigM]=GetOpticsFromSigmaMatrix(TransportedPoints);   % ellypse orientation
        emiGfit(iIter)=max(GetSinglePartEmittance(TransportedPoints,alphaFit(iIter),betaFit(iIter)));  % max 
        MyContours=ExpandMat(MyContours,GenPointsAlongEllypse(alphaFit(iIter),betaFit(iIter),emiGfit(iIter)));
        sigZ(iIter)=sqrt(sigM(1,1)); sigZP(iIter)=sqrt(sigM(2,2));

        % propagated ellypse
        [betaPrp(iIter),alphaPrp(iIter)]=TransportOptics(TMs(:,:,iIter),betaStart,alphaStart);
        emiGprp(iIter)=emiGfit(iIter);
        MyContours=ExpandMat(MyContours,GenPointsAlongEllypse(alphaPrp(iIter),betaPrp(iIter),emiGprp(iIter)));

        % [m,rad] to [mm,mrad]
        TransportedPoints=TransportedPoints*1E3;
        MyContours=MyContours*1E3;

        % calculate histograms
        [nCounts,nCountsX,nCountsY]=Get2dHistograms(TransportedPoints(:,1),TransportedPoints(:,2),xb,yb);
        profiles(:,1+iIter)=nCountsX;
        
        % plot
        if (~ismissing(myTitle))
            Plot2DHistograms(TransportedPoints,nCountsX,nCountsY,xb,yb,"z [mm]","zp [mrad]",MyContours,false);
            sgtitle(sprintf("%s - Transport matrix: %d",myTitle,iIter));
            pause(0.5);
        end
    end
    phisFit=EvaluatePhaseAdvance(TMs,betaStart,alphaStart,betaFit,alphaFit);
    phisPrp=EvaluatePhaseAdvance(TMs,betaStart,alphaStart,betaPrp,alphaPrp);
end
