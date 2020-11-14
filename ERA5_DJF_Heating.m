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

load('ERA5_dsdp_tropics_mean.mat','presplot');
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
lat=lat(latUse);

%file directories
fileDir = '../Raw/ERA5/';
plotDir='../Plots/ERA5/';
qDir = '../Raw/ERA5/q1/';
load('rmm.mat');



%load('ERA5_pcs.mat');
%load the heating profile for the first two EOFs using the time-mean dry static stability

years = 1979:2018;
name = ['../Raw/ERA5/q1/Q1_wtg_' num2str(years(1)) '.mat'];
load(name);
Q1(:,:,:,[1 end])=[];
Q1_tot=zeros(16410,size(Q1,2),size(Q1,3),size(Q1,4));
countst=1;
countnd=size(Q1,1);
Q1_tot(countst:countnd,:,:,:)=Q1;
for i = 2:40
	  name = ['../Raw/ERA5/q1/Q1_wtg_' num2str(years(i)) '.mat'];
clear Q1;
load(name);
Q1(:,:,:,[1 end])=[];
countst=countnd+1;
countnd=countnd+size(Q1,1);
Q1_tot(countst:countnd,:,:,:)=Q1;
disp(years(i));
end
disp('done loading');
%O2shaped=-O2shaped;
%load('ERA5_PCs.mat');

[ntim,np, nlat, nlon] = size(Q1_tot) ;
%np=size(Q1,2);
XX = ones(ntim, 1) ;
%remove the first three harmonics of the yearly cycle
for i = 1:3 ; % Remove the first three harmonics
    XX = [XX sin(2*pi*i*[1:ntim]'/365) cos(2*pi*i*[1:ntim]'/365)] ;
end
disp('remove mean');
XXinv = inv(XX'*XX) ;
latish=-22.5:2.5:22.5;
for i = 1:nlat ;
    for j = 1:nlon ;
        for k = 1:np;
b3= XXinv * ( XX'*squeeze(Q1_tot(:,k,i,j)));
Q1_tot(:,k,i,j) = squeeze(Q1_tot(:,k,i,j))-XX*b3;
        end
    end
	disp(latish(i))
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
disp('get only rmm');
Q1 = Q1_tot(totindex,:,2:end-1,2:end-1);
Q1_rmm=zeros(8,size(Q1,2),size(Q1,3),size(Q1,4));
  for i = 1:8
	    rmmindex=find(rmmphase==i);
Q1_rmm(i,:,:,:)=mean(Q1(rmmindex,:,:,:),1);
  end
  
%%
%calculate the mean top-heaviness angle
%%
%lon=lon(2:end-1);
lon=[lon(73:end); lon(1:72)+360];
Q1_rmm=cat(4,Q1_rmm(:,:,:,73:end),Q1_rmm(:,:,:,1:72))*(86400/1005);
%%
%save the heating rates
save([qDir 'ERA5_q1_norm_RMM_DJF.mat'],'Q1_rmm','lon','lat','presplot');
fout='ERA5_Q1_prime_norm_DJF_MJO.nc';

%create netcdf and add dimensions
nc_create_empty(fout);
nc_adddim(fout,'Phase',0);
varstruct.Name = 'Phases';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'Phase' };
nc_addvar(fout,varstruct);
nc_varput(fout,'Phases',1:8);

nc_adddim(fout,'Pressure',np);
varstruct.Name = 'Pressure';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'Pressure' };
nc_addvar(fout,varstruct);
nc_varput(fout,'Pressure',presplot);


nc_adddim(fout,'Latitude',nlat);
varstruct.Name = 'Latitude';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'Latitude' };
nc_addvar(fout,varstruct);
nc_varput(fout,'Latitude',lat);

nc_adddim(fout,'Longitude',nlon);
varstruct.Name = 'Longitude';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'Longitude' };
nc_addvar(fout,varstruct);
nc_varput(fout,'Longitude',lon);

clear varsturct
%add heating profiles
varstruct.Name = 'Q1_prime';
varstruct.Nctype = 'double';
varstruct.Dimension= {'Phase','Pressure','Latitude','Longitude'};
nc_addvar(fout,varstruct);
nc_varput(fout,varstruct.Name,Q1_rmm);
nc_attput(fout,varstruct.Name','description','apparent heating anomaly due to the first EOF of vertical motion');
nc_attput(fout,varstruct.Name','units','K/day');
nc_dump(fout);
%%



plotting=1;
if plotting
%plot Q1' due to O1
figure('units','inches','Position',[0 0 8 16]),
  title('Q1 anomaly');colormap(redblue);
subplot(8,1,1),
  contourf(lon,lat,squeeze(mean(Q1_rmm(1,:,:,:),2)));
xlabel('phase 1');colorbar;
hold on;
contour(x,y,squeeze(land),'k');
set(gca,'ylim',[-20 20]);
caxis([-2 2])
subplot(8,1,2),
  contourf(lon,lat,squeeze(mean(Q1_rmm(2,:,:,:),2)));
xlabel('phase 2');colorbar;
hold on;
contour(x,y,squeeze(land),'k');
set(gca,'ylim',[-20 20]);
caxis([-2 2])

subplot(8,1,3),
  contourf(lon,lat,squeeze(mean(Q1_rmm(3,:,:,:),2)));
xlabel('phase 3');colorbar;
hold on;
contour(x,y,squeeze(land),'k');
set(gca,'ylim',[-20 20]);
caxis([-2 2])

subplot(8,1,4),
  contourf(lon,lat,squeeze(mean(Q1_rmm(4,:,:,:),2)));
xlabel('phase 4');colorbar;
hold on;
contour(x,y,squeeze(land),'k');
set(gca,'ylim',[-20 20]);
caxis([-2 2])

subplot(8,1,5),
  contourf(lon,lat,squeeze(mean(Q1_rmm(5,:,:,:),2)));
xlabel('phase 5');colorbar;
hold on;
contour(x,y,squeeze(land),'k');
set(gca,'ylim',[-20 20]);
caxis([-2 2])

subplot(8,1,6),
  contourf(lon,lat,squeeze(mean(Q1_rmm(6,:,:,:),2)));
xlabel('phase 6');colorbar;
hold on;
contour(x,y,squeeze(land),'k');
set(gca,'ylim',[-20 20]);
caxis([-2 2])

subplot(8,1,7),
  contourf(lon,lat,squeeze(mean(Q1_rmm(7,:,:,:),2)));
xlabel('phase 7');colorbar;
hold on;
contour(x,y,squeeze(land),'k');
set(gca,'ylim',[-20 20]);
caxis([-2 2])

subplot(8,1,8),
  contourf(lon,lat,squeeze(mean(Q1_rmm(8,:,:,:),2)));
xlabel('phase 8');colorbar;
hold on;
contour(x,y,squeeze(land),'k');
set(gca,'ylim',[-20 20]);
caxis([-2 2])
print(gcf,'-djpeg','ERA5_DJF_rmm_tot_heating.jpg');
 
end
%end
