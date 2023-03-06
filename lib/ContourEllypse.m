function MyEllypse=ContourEllypse(alpha,beta,emiG,dAng)
    if (~exist("dAng","var")), dAng=1; end
    angs=0:dAng:360;
    nEllypses=length(emiG);
    if (length(alpha)==1 && nEllypses>1), alpha=ones(nEllypses,1)*alpha; end
    if (length( beta)==1 && nEllypses>1), beta =ones(nEllypses,1)* beta; end
    MyEllypse=NaN(length(angs),2,nEllypses);
    for ii=1:nEllypses
        MyEllypse(:,1,ii)=cosd(angs);
        MyEllypse(:,2,ii)=sind(angs);
        MyEllypse(:,:,ii)=Norm2Phys(MyEllypse(:,:,ii),beta(ii),alpha(ii),emiG(ii));
    end
end
