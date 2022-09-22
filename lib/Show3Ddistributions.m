function Show3Ddistributions(profiles,addIndex,addLabel,myCaseNames,myTitle)
    nIter=size(profiles,3);
    [nRows,nCols]=GetNrowsNcols(nIter);
    ff=figure();
    ff.Position(1:2)=[1000 100]; % figure at lower-left corner of screen
    ff.Position(3)=ff.Position(3)*nCols; % larger figure
    ff.Position(4)=ff.Position(4)*nRows; % larger figure
    t = tiledlayout(nRows,nCols,'TileSpacing','Compact','Padding','Compact'); % minimise whitespace around plots
    for iIter=1:nIter
        nexttile; % subplot(nRows,nCols,iIter);
        PlotSpectra(profiles(:,:,iIter),false,addIndex,addLabel);
        title(myCaseNames(iIter));
        xlabel("position [mm]");
        zlabel("Counts []");
    end
    sgtitle(myTitle);
end
