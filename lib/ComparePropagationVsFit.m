function ComparePropagationVsFit(Xs,xLab,alphaFit,betaFit,emiGfit,phisFit,alphaPrp,betaPrp,emiGprp,phisPrp,sigZ,sigZP,myTitles)
    for iFig=1:size(betaFit,2)
        ff=figure();
        ff.Position(1:2)=[1000 50]; % figure at lower-left corner of screen
        ff.Position(3:4)=ff.Position(3:4)*3; % larger figure
        t = tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact'); % minimise whitespace around plots
        ii=0;
        % - beta
        ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
        plot(Xs,betaFit(:,iFig),"o-",Xs,betaPrp(:,iFig),"*-");
        grid(); xlabel(xLab); ylabel("\beta [m]"); legend("Fit","Propagated","Location","best");
        % - alpha
        ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
        plot(Xs,alphaFit(:,iFig),"o-",Xs,alphaPrp(:,iFig),"*-");
        grid(); xlabel(xLab); ylabel("\alpha []"); legend("Fit","Propagated","Location","best");
        % - emiG
        ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
        plot(Xs,emiGfit(:,iFig)*1E6,"o-",Xs,emiGprp(:,iFig)*1E6,"*-");
        grid(); xlabel(xLab); ylabel("\epsilon [\mum]"); legend("Fit","Propagated","Location","best");
        % - sig_z
        ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
        plot(Xs,sigZ(:,iFig)*1E3,"o-",Xs,sqrt(betaFit(:,iFig).*emiGfit(:,iFig))*1E3,"*-");
        grid(); xlabel(xLab); ylabel("\sigma_{z} [mm]"); legend("Stats","Fit","Location","best");
        % - sig_zp
        ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
        plot(Xs,sigZP(:,iFig)*1E3,"o-",Xs,sqrt((1+alphaFit(:,iFig).^2)./betaFit(:,iFig).*emiGfit(:,iFig))*1E3,"*-");
        grid(); xlabel(xLab); ylabel("\sigma_{zp} [mrad]"); legend("Stats","Fit","Location","best");
        % - phi
        ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
        plot(Xs,phisFit(:,iFig)/pi*180,"o-",Xs,phisPrp(:,iFig)/pi*180,"*-");
        grid(); xlabel(xLab); ylabel("\phi [deg]"); legend("Fit","Propagated","Location","best");
        % - ratio sig_z
        ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
        plot(Xs,sigZ(:,iFig)./sqrt(betaFit(:,iFig).*emiGfit(:,iFig)),"o-");
        grid(); xlabel(xLab); ylabel("R_{\sigma_{z}} [mm]");
        % - ratio sig_zp
        ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
        plot(Xs,sigZP(:,iFig)./sqrt((1+alphaFit(:,iFig).^2)./betaFit(:,iFig).*emiGfit(:,iFig)),"o-");
        grid(); xlabel(xLab); ylabel("R_{\sigma_{z}} [mm]");
        % - ratio phis
        ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
        plot(Xs,phisFit(:,iFig)./phisPrp(:,iFig),"o-");
        grid(); xlabel(xLab); ylabel("R_{\sigma_{z}} [mm]");
        %
        linkaxes(ax,"x");
        sgtitle(myTitles(iFig));
    end
end
