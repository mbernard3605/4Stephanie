%This script calculates the Q1' heating from the first two PCs of vertical motion over DJF for the 40
%year climatology from ERA5.  

%Add the paths for the netcdf scripts 
   addpath('~/matlab_scripts/mexcdf/mexnc');
   addpath('~/matlab_scripts/mexcdf/snctools');
%load objects needed for plotting
load('anglecolormap2.mat');
testmap=anglemap;
testmap=[testmap(33:end,:); testmap(1:32,:)];
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
lats = nc_varget(landFile,'latitude');

latUse=find(lats>=-22.5 & lats<=22.5);
lat=lats(latUse);
lon=[lon(72:end); lon(1:71)+360];
[x,y]=meshgrid(lon,lat);

%file directories
fileDir = '../Raw/ERA5/';
plotDir='../Plots/ERA5/';
load('rmm.mat');


%load('ERA5_pcs.mat');
%load the heating profile for the first two EOFs using the time-mean dry static stability

%O2shaped=-O2shaped;
load('ERA5_PCs.mat');
disp('remove mean');

O1tot(:,:,[1 end]) = [];
O2tot(:,:,[1 end]) = [];


[ntim, nlat, nlon] = size(O1tot) ;
XX = ones(ntim, 1) ;
%remove the first three harmonics of the yearly cycle
%O1full_detrend=zeros(size(O1tot));O2full_detrend=O1full_detrend;
for i = 1:3 ; % Remove the first three harmonics
    XX = [XX sin(2*pi*i*[1:ntim]'/365) cos(2*pi*i*[1:ntim]'/365)] ;
end

XXinv = inv(XX'*XX) ;
for i = 1:nlat ;
    for j = 1:nlon ;
            b1 = XXinv * (XX'*squeeze(O1tot(:,i,j))) ;
            b2 = XXinv * (XX'*squeeze(O2tot(:,i,j))) ;
            O1tot(:,i,j) = squeeze(O1tot(:,i,j)) - XX*b1 ;
            O2tot(:,i,j) = squeeze(O2tot(:,i,j)) - XX*b2 ;
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
%calculate the heating rates at each location and each time

%calculate the mean top-heaviness angle
%%
[nT,nLat,nLon]=size(O1tot);
O1_prime_rmm=zeros(8,nLat,nLon);
O2_prime_rmm=zeros(8,nLat,nLon);
disp('compositing');			  
%Composite by MJO index
for i = 1:8
	  rmmindex=find(rmmphase==i);
O1_prime_rmm(i,:,:)=mean(O1tot(rmmindex,:,:),1);
O2_prime_rmm(i,:,:)=mean(O2tot(rmmindex,:,:),1);
end
O1_prime_rmm=cat(3,O1_prime_rmm(:,:,72:end),O1_prime_rmm(:,:,1:71));
O2_prime_rmm=cat(3,O2_prime_rmm(:,:,72:end),O2_prime_rmm(:,:,1:71));
angle_prime_rmm=atan2d(O2_prime_rmm,O1_prime_rmm);
%%
disp('saving');
%save the heating rates
save(['ERA5_angle_decom_RMM_DJF.mat'],'angle_prime_rmm','O1_prime_rmm','O2_prime_rmm','lon','lat');
%%
plotting=1;
if plotting
%plot Q1' due to O1
figure('units','inches','Position',[0 0 8 16]),
colormap(redblue);
  title('Q1 anomaly due to EOF 1');
subplot(8,1,1),
  contourf(lon,lat,squeeze(O1_prime_rmm(1,:,:)));
xlabel('phase 1');colorbar;
hold on;
contour(lon,lats,squeeze(land));
subplot(8,1,2),
  contourf(lon,lat,squeeze(O1_prime_rmm(2,:,:)));
xlabel('phase 2');colorbar;
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,3),
  contourf(lon,lat,squeeze(O1_prime_rmm(3,:,:)));
xlabel('phase 3');colorbar;
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,4),
  contourf(lon,lat,squeeze(O1_prime_rmm(4,:,:)));
xlabel('phase 4');colorbar;
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,5),
  contourf(lon,lat,squeeze(O1_prime_rmm(5,:,:)));
xlabel('phase 5');colorbar;
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,6),
  contourf(lon,lat,squeeze(O1_prime_rmm(6,:,:)));
xlabel('phase 6');colorbar;
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,7),
  contourf(lon,lat,squeeze(O1_prime_rmm(7,:,:)));
xlabel('phase 7');colorbar;
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,8),
  contourf(lon,lat,squeeze(O1_prime_rmm(8,:,:)));
xlabel('phase 8');colorbar;
hold on;
contour(lon,lats,squeeze(land));


print(gcf,'-djpeg',[plotDir 'DJF_RMM_o1.jpg']);


figure('units','inches','Position',[0 0 8 16]),
  colormap(redblue);
  title('Q1 anomaly due to EOF 2');
subplot(8,1,1),
  contourf(lon,lat,squeeze(O2_prime_rmm(1,:,:)));
xlabel('phase 1');colorbar;
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,2),
  contourf(lon,lat,squeeze(O2_prime_rmm(2,:,:)));
xlabel('phase 2');colorbar;
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,3),
  contourf(lon,lat,squeeze(O2_prime_rmm(3,:,:)));
xlabel('phase 3');colorbar;
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,4),
  contourf(lon,lat,squeeze(O2_prime_rmm(4,:,:)));
xlabel('phase 4');colorbar;
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,5),
  contourf(lon,lat,squeeze(O2_prime_rmm(5,:,:)));
xlabel('phase 5');colorbar;
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,6),
  contourf(lon,lat,squeeze(O2_prime_rmm(6,:,:)));
xlabel('phase 6');colorbar;
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,7),
  contourf(lon,lat,squeeze(O2_prime_rmm(7,:,:)));
xlabel('phase 7');colorbar;
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,8),
  contourf(lon,lat,squeeze(O2_prime_rmm(8,:,:)));
xlabel('phase 8');colorbar;
hold on;
contour(lon,lats,squeeze(land));


print(gcf,'-djpeg',[plotDir 'DJF_RMM_o2.jpg']);



figure('units','inches','Position',[0 0 8 16]),
  colormap(testmap);
  title('Q1 anomaly due to EOF 1+2');
subplot(8,1,1),
  contourf(lon,lat,squeeze(angle_prime_rmm(1,:,:)));
xlabel('phase 1');colorbar;
caxis([-180 180]);
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,2),
  contourf(lon,lat,squeeze(angle_prime_rmm(2,:,:)));
xlabel('phase 2');colorbar;
caxis([-180 180]);
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,3), 
  contourf(lon,lat,squeeze(angle_prime_rmm(3,:,:)));
xlabel('phase 3');colorbar;
caxis([-180 180]);
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,4),
  contourf(lon,lat,squeeze(angle_prime_rmm(4,:,:)));
xlabel('phase 4');colorbar;
caxis([-180 180]);
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,5),
  contourf(lon,lat,squeeze(angle_prime_rmm(5,:,:)));
xlabel('phase 5');colorbar;
caxis([-180 180]);
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,6),
  contourf(lon,lat,squeeze(angle_prime_rmm(6,:,:)));
xlabel('phase 6');colorbar;
caxis([-180 180]);
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,7),
  contourf(lon,lat,squeeze(angle_prime_rmm(7,:,:)));
xlabel('phase 7');colorbar;
caxis([-180 180]);
hold on;
contour(lon,lats,squeeze(land));

subplot(8,1,8),
  contourf(lon,lat,squeeze(angle_prime_rmm(8,:,:)));
xlabel('phase 8');colorbar;
caxis([-180 180]);
hold on;
contour(lon,lats,squeeze(land));


print(gcf,'-djpeg',[plotDir 'DJF_RMM_angle.jpg']);


 
end
%end
