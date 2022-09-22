function MyEllypse=SampleEllypse(alpha,beta,emiG,nEl)
    if (~exist("nEl","var")), nEl=360; end
    nEllypses=size(emiG,2);
    MyEllypse=NaN(nEl+1,2,nEllypses);
    for ii=1:nEllypses
        MyEllypse(1:end-1,1,ii)=cosd(linspace(0,360,nEl));
        MyEllypse(1:end-1,2,ii)=sind(linspace(0,360,nEl));
        MyEllypse(end,:,ii)=MyEllypse(1,:,ii);
        [MyEllypse(:,:,ii)]=Norm2Phys(MyEllypse(:,:,ii),beta(ii),alpha(ii),emiG(ii));
    end
end
