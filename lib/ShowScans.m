function ShowScans(Xs,xLab,alpha,beta,emiG,phis,sigZ,sigZP,myLegend,myTitle)
    ff=figure();
    ff.Position(1:2)=[1000 50]; % figure at lower-left corner of screen
    ff.Position(3:4)=ff.Position(3:4)*3; % larger figure
    nSeries=size(beta,2);
    t = tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact'); % minimise whitespace around plots
    if ( nSeries==1 )
        cm=zeros(1,3);
    else
        cm=colormap(parula(nSeries));
    end
    ii=0;
    % - beta
    ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
    for iSeries=1:nSeries
        if ( iSeries>1 ), hold on; end
        plot(Xs,beta(:,iSeries),".-","Color",cm(iSeries,:));
        grid(); xlabel(xLab); ylabel("\beta [m]");
    end
    % - alpha
    ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
    for iSeries=1:nSeries
        if ( iSeries>1 ), hold on; end
        plot(Xs,alpha(:,iSeries),".-","Color",cm(iSeries,:));
        grid(); xlabel(xLab); ylabel("\alpha []");
    end
    % - emiG
    ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
    for iSeries=1:nSeries
        if ( iSeries>1 ), hold on; end
        plot(Xs,emiG(:,iSeries)*1E6,".-","Color",cm(iSeries,:));
        grid(); xlabel(xLab); ylabel("\epsilon [\mum]");
    end
    % - sig_z
    ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
    for iSeries=1:nSeries
        if ( iSeries>1 ), hold on; end
        plot(Xs,sigZ(:,iSeries)*1E3,".-","Color",cm(iSeries,:));
        grid(); xlabel(xLab); ylabel("\sigma_{z} [mm]");
    end
    % - sig_zp
    ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
    for iSeries=1:nSeries
        if ( iSeries>1 ), hold on; end
        plot(Xs,sigZP(:,iSeries)*1E3,".-","Color",cm(iSeries,:));
        grid(); xlabel(xLab); ylabel("\sigma_{zp} [mrad]");
    end
    % - phi
    ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
    for iSeries=1:nSeries
        if ( iSeries>1 ), hold on; end
        plot(Xs,phis(:,iSeries)/pi*180,".-","Color",cm(iSeries,:));
        grid(); xlabel(xLab); ylabel("\phi [deg]");
    end
    % - legend
    ii=ii+1; ax(ii)=nexttile; % subplot(3,3,ii);
    for iSeries=1:nSeries
        if ( iSeries>1 ), hold on; end
        plot(NaN(),NaN(),".-","Color",cm(iSeries,:));
    end
    legend(myLegend,"Location","best");
    % - general
    linkaxes(ax,"x");
    sgtitle(myTitle);
end
