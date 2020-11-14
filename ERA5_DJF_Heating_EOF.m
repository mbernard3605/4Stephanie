%This script calculates the Q1' heating from the first two PCs of vertical motion over DJF for the 40
%year climatology from ERA5.  

%Add the paths for the netcdf scripts 
   addpath('~/matlab_scripts/mexcdf/mexnc');
   addpath('~/matlab_scripts/mexcdf/snctools');
%load objects needed for plotting
load('anglecolormap2.mat');
testmap=anglemap;
testmap=[testmap(33:end,:) testmap(1:32,:)];
load('coasts.mat')
Cp=1005;
g=9.81;
Lv=2.5e6;


%base directory 
baseDir='../Raw/ERA5/';

%needed if land-masking otherwise is needed to plot contours
 landFile=[baseDir 'ERA5_land.nc'];
 land = nc_varget(landFile,'lsm');
% sea=~land;
% land(find(landMask))=NaN;

%land = zeros(size(land));
lon = nc_varget(landFile,'longitude');
lat = nc_varget(landFile,'latitude');
[x,y]=meshgrid([lon(73:end);lon(1:72)+360],lat);                                                        
land=cat(3,land(:,:,73:end),land(:,:,1:72));

latUse=find(lat>=-20 & lat<=20);
lats=lat(latUse);
lons=[lon(72:end); lon(1:71)+360];


%file directories
fileDir = '../Raw/ERA5/';
plotDir='../Plots/ERA5/';
qDir = '../Raw/ERA5/q1/';
load('rmm.mat');
load([qDir 'ERA5_q1_tot_mean.mat']);


%load('ERA5_pcs.mat');
%load the heating profile for the first two EOFs using the time-mean dry static stability

%O2shaped=-O2shaped;
load('ERA5_PCs.mat');
disp('remove mean');
O1tot(:,:,[1 end]) = [];
O2tot(:,:,[1 end]) = [];
Q1_O1_tot(:,:,[1 end]) = [];
Q1_O2_tot(:,:,[1 end]) = [];
lon=lon(2:end-1); 
[ntim, nlat, nlon] = size(O1tot) ;
XX = ones(ntim, 1) ;
%remove the first three harmonics of the yearly cycle
O1full_detrend=zeros(size(O1tot));O2full_detrend=O1full_detrend;
for i = 1:3 ; % Remove the first three harmonics
    XX = [XX sin(2*pi*i*[1:ntim]'/365) cos(2*pi*i*[1:ntim]'/365)] ;
end

XXinv = inv(XX'*XX) ;
for i = 1:nlat ;
    for j = 1:nlon ;
            b1 = XXinv * (XX'*squeeze(O1tot(:,i,j))) ;
            b2 = XXinv * (XX'*squeeze(O2tot(:,i,j))) ;
            O1full_detrend(:,i,j) = squeeze(O1tot(:,i,j)) - XX*b1 ;
            O2full_detrend(:,i,j) = squeeze(O2tot(:,i,j)) - XX*b2 ;
    end
end



%months that we are selecting
    months=[12, 1, 2];

%get the subset of time that has |rmm|^2>=1
starttime=1979;
endtime=2019;
startindex=find(rmmsave(:,1)==starttime);
endindex=find(rmmsave(:,1)==endtime);
rmm=rmmsave(startindex(1):endindex(1)-1,:);
RMMindex=find(rmm(:,5)>=1);
timeindex=find(rmm(:,2)==months(1) | rmm(:,2)==months(2) | rmm(:,2)==months(3) );
totindex=intersect(RMMindex,timeindex);
time=rmm(totindex,1:3);
rmmphase=rmm(totindex,4);



%%
    Q1_O1_prime=zeros(size(totindex,1),size(Q1_O1_tot,1),size(Q1_O1_tot,2),size(Q1_O1_tot,3)-2); 
    Q1_O2_prime=zeros(size(totindex,1),size(Q1_O2_tot,1),size(Q1_O2_tot,2),size(Q1_O2_tot,3)-2);
%calculate the heating rates at each location and each time
O1full_detrend_big=permute(repmat(O1full_detrend(totindex,:,:),[1,1,1,size(Q1_O1_tot,1)]),[1 4 2 3]);
O2full_detrend_big=permute(repmat(O2full_detrend(totindex,:,:),[1,1,1,size(Q1_O2_tot,1)]),[1 4 2 3]);
Q1_O1_big=permute(repmat(Q1_O1_tot,[1,1,1,length(totindex)]),[4 1 2 3]);
Q1_O2_big=permute(repmat(Q1_O2_tot,[1,1,1,length(totindex),]),[4 1 2 3]);

Q1_O1_prime=O1full_detrend_big.*Q1_O1_big/(Cp/86400);
Q1_O2_prime=O2full_detrend_big.*Q1_O2_big/(Cp/86400);

clear O1full_detrend_big O2full_detrend_big Q1_O1_big Q1_O2_big
%calculate the mean top-heaviness angle
%  ERA5angle=atan2d(squeeze(nanmean(O2tot,1)),squeeze(nanmean(O1tot)));
%%

[nT,nP,nLat,nLon]=size(Q1_O1_prime);
Q1_O1_prime_rmm=zeros(8,nP,nLat,nLon);
Q1_O2_prime_rmm=zeros(8,nP,nLat,nLon);
disp('compositing');			  
%Composite by MJO index
for i = 1:8
	  rmmindex=find(rmmphase==i);
Q1_O1_prime_rmm(i,:,:,:)=mean(Q1_O1_prime(rmmindex,:,:,:),1);
Q1_O2_prime_rmm(i,:,:,:)=mean(Q1_O2_prime(rmmindex,:,:,:),1);
				  
end
Q1_O1_prime_rmm=cat(4,Q1_O1_prime_rmm(:,:,:,72:end),Q1_O1_prime_rmm(:,:,:,1:71));
Q1_O2_prime_rmm=cat(4,Q1_O2_prime_rmm(:,:,:,72:end),Q1_O2_prime_rmm(:,:,:,1:71));
Q1_tot_prime_rmm=Q1_O1_prime_rmm+Q1_O2_prime_rmm;

%%
disp('saving');
%save the heating rates
save([qDir 'ERA5_q1_decom_RMM_DJF.mat'],'Q1_O1_prime_rmm','Q1_O2_prime_rmm','lon','lat','presplot');
fout='ERA5_Q1_prime_O1_O2_DJF_MJO.nc';

%create netcdf and add dimensions
nc_create_empty(fout);
nc_adddim(fout,'Phase',0);
varstruct.Name = 'Phases';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'Phase' };
nc_addvar(fout,varstruct);
nc_varput(fout,'Phases',1:8);

nc_adddim(fout,'Pressure',nP);
varstruct.Name = 'Pressure';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'Pressure' };
nc_addvar(fout,varstruct);
nc_varput(fout,'Pressure',presplot);


nc_adddim(fout,'Latitude',nLat);
varstruct.Name = 'Latitude';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'Latitude' };
nc_addvar(fout,varstruct);
nc_varput(fout,'Latitude',lat);

nc_adddim(fout,'Longitude',nLon);
varstruct.Name = 'Longitude';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'Longitude' };
nc_addvar(fout,varstruct);
nc_varput(fout,'Longitude',lons);

clear varsturct
%add heating profiles
varstruct.Name = 'Q1_prime_O1';
varstruct.Nctype = 'double';
varstruct.Dimension= {'Phase','Pressure','Latitude','Longitude'};
nc_addvar(fout,varstruct);
nc_varput(fout,varstruct.Name,Q1_O1_prime_rmm);
nc_attput(fout,varstruct.Name','description','apparent heating anomaly due to the first EOF of vertical motion');
nc_attput(fout,varstruct.Name','units','K/day');

varstruct.Name = 'Q1_prime_O2';
varstruct.Nctype = 'double';
varstruct.Dimension = {'Phase','Pressure','Latitude','Longitude'};
nc_addvar(fout,varstruct);
nc_varput(fout,varstruct.Name,Q1_O2_prime_rmm);
nc_attput(fout,varstruct.Name','description','apparent heating anomaly due to the second EOF of vertical motion');
nc_attput(fout,varstruct.Name','units','K/day');

varstruct.Name = 'Q1_prime_O1_O2';
varstruct.Nctype = 'double';
varstruct.Dimension={'Phase','Pressure','Latitude','Longitude'};
nc_addvar(fout,varstruct);
nc_varput(fout,varstruct.Name,Q1_tot_prime_rmm);
nc_attput(fout,varstruct.Name','description','apparent heating anomaly due to the first and second EOFs of vertical motion');
nc_attput(fout,varstruct.Name','units','K/day');
disp('plotting');
nc_dump(fout);
%%
plotting=1;
if plotting
%plot Q1' due to O1
figure('units','inches','Position',[0 0 8 16]),
colormap(redblue);
  title('Q1 anomaly due to EOF 1');
subplot(8,1,1),
  contourf(lons,lat,squeeze(mean(Q1_O1_prime_rmm(1,:,:,:),2)));
xlabel('phase 1');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,2),
  contourf(lons,lat,squeeze(mean(Q1_O1_prime_rmm(2,:,:,:),2)));
xlabel('phase 2');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,3),
  contourf(lons,lat,squeeze(mean(Q1_O1_prime_rmm(3,:,:,:),2)));
xlabel('phase 3');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,4),
  contourf(lons,lat,squeeze(mean(Q1_O1_prime_rmm(4,:,:,:),2)));
xlabel('phase 4');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,5),
  contourf(lons,lat,squeeze(mean(Q1_O1_prime_rmm(5,:,:,:),2)));
xlabel('phase 5');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,6),
  contourf(lons,lat,squeeze(mean(Q1_O1_prime_rmm(6,:,:,:),2)));
xlabel('phase 6');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,7),
  contourf(lons,lat,squeeze(mean(Q1_O1_prime_rmm(7,:,:,:),2)));
xlabel('phase 7');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,8),
  contourf(lons,lat,squeeze(mean(Q1_O1_prime_rmm(8,:,:,:),2)));
xlabel('phase 8');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);


print(gcf,'-djpeg',[plotDir 'DJF_RMM_heating_o1.jpg']);


figure('units','inches','Position',[0 0 8 16]),
colormap(redblue);

  title('Q1 anomaly due to EOF 2');
subplot(8,1,1),
  contourf(lons,lat,squeeze(mean(Q1_O2_prime_rmm(1,:,:,:),2)));
xlabel('phase 1');colorbar;
caxis([-.5 .5]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,2),
  contourf(lons,lat,squeeze(mean(Q1_O2_prime_rmm(2,:,:,:),2)));
xlabel('phase 2');colorbar;
caxis([-.5 .5]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,3),
  contourf(lons,lat,squeeze(mean(Q1_O2_prime_rmm(3,:,:,:),2)));
xlabel('phase 3');colorbar;
caxis([-.5 .5]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,4),
  contourf(lons,lat,squeeze(mean(Q1_O2_prime_rmm(4,:,:,:),2)));
xlabel('phase 4');colorbar;
caxis([-.5 .5]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,5),
  contourf(lons,lat,squeeze(mean(Q1_O2_prime_rmm(5,:,:,:),2)));
xlabel('phase 5');colorbar;
caxis([-.5 .5]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,6),
  contourf(lons,lat,squeeze(mean(Q1_O2_prime_rmm(6,:,:,:),2)));
xlabel('phase 6');colorbar;
caxis([-.5 .5]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,7),
  contourf(lons,lat,squeeze(mean(Q1_O2_prime_rmm(7,:,:,:),2)));
xlabel('phase 7');colorbar;
caxis([-.5 .5]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,8),
  contourf(lons,lat,squeeze(mean(Q1_O2_prime_rmm(8,:,:,:),2)));
xlabel('phase 8');colorbar;
caxis([-.5 .5]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);


print(gcf,'-djpeg',[plotDir 'DJF_RMM_heating_o2.jpg']);



figure('units','inches','Position',[0 0 8 16]),
  colormap(redblue);
  title('Q1 anomaly due to EOF 1+2');
subplot(8,1,1),
  contourf(lons,lat,squeeze(mean(Q1_tot_prime_rmm(1,:,:,:),2)));
xlabel('phase 1');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,2),
  contourf(lons,lat,squeeze(mean(Q1_tot_prime_rmm(2,:,:,:),2)));
xlabel('phase 2');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,3), 
  contourf(lons,lat,squeeze(mean(Q1_tot_prime_rmm(3,:,:,:),2)));
xlabel('phase 3');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,4),
  contourf(lons,lat,squeeze(mean(Q1_tot_prime_rmm(4,:,:,:),2)));
xlabel('phase 4');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,5),
  contourf(lons,lat,squeeze(mean(Q1_tot_prime_rmm(5,:,:,:),2)));
xlabel('phase 5');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,6),
  contourf(lons,lat,squeeze(mean(Q1_tot_prime_rmm(6,:,:,:),2)));
xlabel('phase 6');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,7),
  contourf(lons,lat,squeeze(mean(Q1_tot_prime_rmm(7,:,:,:),2)));
xlabel('phase 7');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);

subplot(8,1,8),
  contourf(lons,lat,squeeze(mean(Q1_tot_prime_rmm(8,:,:,:),2)));
xlabel('phase 8');colorbar;
caxis([-2 2]);
hold on;                                                                                                
contour(x,y,squeeze(land),'k');                                                                         
set(gca,'ylim',[-20 20]);


print(gcf,'-djpeg',[plotDir 'DJF_RMM_heating_o1_o2.jpg']);


 
end
%end
