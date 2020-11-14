addpath('~/matlab_scripts/mexcdf/mexnc');
addpath('~/matlab_scripts/mexcdf/snctools');

baseDir='../Raw/ERA5/';
inDir= [baseDir 'w/'];
plotDir='../Plots/';

dataDir=dir([inDir '*.nc']);


landFile=[baseDir 'ERA5_land.nc'];
lonland = nc_varget(landFile,'longitude');
latland = nc_varget(landFile,'latitude');
 latUse=find(latland<=20 & latland>=-20);

land = nc_varget(landFile,'lsm');
land=land(:,latUse,:);
%land = zeros(size(land));
% latland=latland(latUse);

lsm=find(~reshape(logical(land),1,[]));
landsize=numel(lsm);

years=1979:1:2018;
months=1:1:12;

dataFull=zeros(14610,28,17,144);
fillIndex=1;

tic
for yearIndex = 1:length(years)
    for monthIndex = 1:length(months)
        fileName=sprintf('ERA5_w_%d_%d.nc',years(yearIndex),months(monthIndex));
        data=nc_varget([inDir fileName],'w');
        for avgIndex=1:size(data,1)/4
        dataFull(fillIndex,:,:,:)=mean(data(1+4*(avgIndex-1):4*avgIndex,:,latUse,:),1);
        fillIndex=fillIndex+1;
        end
    end
    disp(yearIndex)
end
toc

lon=nc_varget([inDir fileName],'longitude');
lat=nc_varget([inDir fileName],'latitude');
lat=lat(latUse);
level=nc_varget([inDir fileName],'level')*100;

ilat=length(lat);
ilon=length(lon);
itim=fillIndex-1;

dp=diff(level);
weights=dp./sum(dp);
presEnd=find(level==10000);
level=level(presEnd:end);
ilvl=length(level);
[x,y]=meshgrid(lon,lat);

dataFull=dataFull(1:fillIndex-1,presEnd:end,:,:);



dataFull=reshape(permute(dataFull,[2 1 3 4]),ilvl,[]);

dataShape=zeros(ilvl,itim*landsize);
for i = 1:itim
   dataShape(:,1+(i-1)*landsize:i*landsize)=dataFull(:,itim*lsm-itim+i);
end

dataMean=mean(dataShape,2);
dataShape=(dataShape-repmat(dataMean,1,size(dataShape,2))).*repmat(weights.^0.5,1,size(dataShape,2));

Nstar=40;

C=dataShape*dataShape'/(size(dataShape,2)-1);
[lds, lam]=eig(C);
l=diag(lam);
[lam,k]=sort(l');
lam = fliplr(lam);
lds = fliplr(lds);
eof=lds./repmat(weights.^0.5,1,27);


pcs=dataShape'*lds;
[m,n]=size(dataShape);
save('ERA5_EOFS.mat','lds','dataMean','level','weights','lam');
%save('NCEP_PCS.mat','lds');
EOF = lds./repmat(weights.^0.5,1,size(pcs,2));
O1 = pcs(:,1);
O2 = pcs(:,2);
O3 = pcs(:,3);
O4 = pcs(:,4);
O5 = pcs(:,5);
save('ERA5_O1.mat','O1','itim','ilat','ilon','lam');
save('ERA5_O2.mat','O2','itim','ilat','ilon','lam');
save('ERA5_Ox.mat','O3','O4','O5','itim','ilat','ilon','lam');

%Variance explain
%Variance explained plot
% figure,
% errorbar(1:10, per(1:10), sqrt(2/Nstar)*per(1:10)) ;
% set(gca, 'XTick', 1:10) ;
% xlabel('Mode Number') ;
% ylabel('Percent Variance Explained') ;

figure(1),
plot(EOF(:,1:2),level,'linewidth',2);
set(gca,'ydir','reverse','ylim',[10000 100000]);
title('EOFs of vertical motion for ERA5 daily');
legend('EOF1','EOF2');
print(gcf,'-djpeg','-r300',[plotDir 'ERA5_EOF.jpg']);
%%
O1shaped=NaN(size(dataFull,2),1);
O2shaped=O1shaped;

for i = 1:itim
    O1shaped(itim*lsm-itim+i)=O1(1+(i-1)*landsize:i*landsize);
    O2shaped(itim*lsm-itim+i)=O2(1+(i-1)*landsize:i*landsize);
end

O1shaped=reshape(O1shaped,itim,ilat,ilon);
O2shaped=reshape(O2shaped,itim,ilat,ilon);

O1shaped=cat(3,O1shaped(:,:,73:end),O1shaped(:,:,1:72));
O2shaped=cat(3,O2shaped(:,:,73:end),O2shaped(:,:,1:72));

convMask=squeeze(O1shaped)<=0;
convMaskfull=squeeze(mean(O1shaped))<=0;

O1conv=nan(size(O1shaped));O2conv=nan(size(O1shaped));
O1conv(convMask)=O1shaped(convMask);
O2conv(convMask)=O2shaped(convMask);

lon_map=[lon(73:end);lon(1:72)+360];
[xmap,ymap]=meshgrid(lon_map,lat);

figure(2),
contourf(xmap,ymap,squeeze(mean(O1shaped)));
title('Mean O1 ERA5');
colormap(jet);
axis('equal');
print(gcf,'-djpeg','-r300',[plotDir 'ERA5_O1map.jpg']);

figure(3),
contourf(xmap,ymap,squeeze(mean(O2shaped)));
colormap(jet);
title('Mean O2 ERA5');
axis('equal');
print(gcf,'-djpeg','-r300',[plotDir 'ERA5_O2map.jpg']);

figure(4),
contourf(xmap,ymap,squeeze(nanmean(O1conv)));
colormap(jet);
title('Mean O1 (O1<0) ERA5');
axis('equal');
print(gcf,'-djpeg','-r300',[plotDir 'ERA5_O1map_conv.jpg']);

figure(5),
contourf(xmap,ymap,squeeze(nanmean(O2conv)));
colormap(jet);
title('Mean O2 (O1<0) ERA5');
axis('equal');
print(gcf,'-djpeg','-r300',[plotDir 'ERA5_O2map_conv.jpg']);

OSS=squeeze(nanmean(O2conv))./squeeze(nanmean(O1conv))*sqrt(lam(1)/lam(2));
OSS=[OSS(:,73:end) OSS(:,72:end)];

figure(6),
contourf(xmap,ymap,squeeze(nanmean(O2conv))./squeeze(nanmean(O1conv))*sqrt(lam(1)/lam(2)),-1:.05:1);
colormap(redblue);
colorbar;
title('Mean O2/O1 ratio ERA5');
axis('equal');
hold on;
contour(xmap,ymap,squeeze(mean(O1shaped)),[0,0],'k','linewidth',2);
print(gcf,'-djpeg','-r300',[plotDir 'ERA5_O2O1map_conv.jpg']);

figure(7),
contourf(xmap,ymap,squeeze(nanmean(O2shaped))./squeeze(nanmean(O1shaped))*sqrt(lam(1)/lam(2)),-1:0.05:1);
colormap(redblue);
axis('equal');

% 
% save('ERA5_eofs.mat','EOF');
 save('ERA5_pcs.mat','O1shaped','O2shaped')
% save('ERA5_pcsraw.mat','O1','O2');
