%function ERA5_Heating()
  %if running as a compiled program, this is needs to be turned off and these folders need to be added at compliation
if ~isdeployed
addpath('~/iris-home/matlab_scripts/mexcdf/mexnc');
addpath('~/iris-home/matlab_scripts/mexcdf/snctools');
end
%load('NCEP_EOFS.mat','omegaMean');

load('ERA5_EOFS.mat');

%define constants
g=9.81;
earthRad=6.371e6;
Cp=1005;

%directories where all the variables and output are stored
baseDir='../Raw/ERA5/';
tDir = [baseDir 't/'];
zDir = [baseDir 'z/'];
qDir = [baseDir 'q1/'];
plotDir='../Plots/ERA5/';


%load the landmask, this is not currently being used to mask the data
landFile=[baseDir 'ERA5_land.nc'];
landMask = nc_varget(landFile,'lsm');
land=ones(size(landMask));
land(find(landMask))=NaN;

%land = zeros(size(land));
lonland = nc_varget(landFile,'longitude');
latland = nc_varget(landFile,'latitude');

%define the time-span
years=1979:1:2018;
months=1:1:12;
started=0;
%allocate the variables


for yearsIndex = 1:length(years)
    startIndex=1;
    outFile=[qDir 'Q1_eof_' num2str(years(yearsIndex)) '.nc'];
if ~exist(outFile)
    if ~mod(years(yearsIndex),4)
        baseFile='../baseleap_ERA5.nc';
        days=366;
        TempFull=zeros(days,28,19,144);
        GeopFull=TempFull;
        timeFull=zeros(days,1);
    else
        baseFile='../base_ERA5.nc';
        days=365;
        TempFull=zeros(days,28,19,144);
        GeopFull=TempFull;
        timeFull=zeros(days,1);
    end
copyfile(baseFile,outFile);
end

for monthIndex = 1:length(months);

            fileAppend=sprintf('_%d_%d.nc',years(yearsIndex),months(monthIndex));
time = nc_varget([tDir 'ERA5_t' fileAppend],'time');

if ~started
lon=nc_varget([tDir 'ERA5_t' fileAppend],'longitude');
lat=nc_varget([tDir 'ERA5_t' fileAppend],'latitude');
latUse = find(lat<=22.5 & lat>=-22.5); 
lon = [lon(end); lon; lon(1)];
lat=lat(latUse);
level=nc_varget([tDir 'ERA5_t' fileAppend],'level')*100;
presEnd=find(level==10000);
presEnd=1;
dp=diff(level);
weights=dp./sum(dp);

pres=level(presEnd:end);
started=1;

S_tot=zeros(size(dp,2),size(dp,3),size(dp,4));
[xplot,yplot]=meshgrid(lon(2:end-1),lat(2:end-1));
presplot=(pres(2:end)+pres(1:end-1))/2;

end




Temp=(nc_varget([tDir 'ERA5_t' fileAppend],'t'));
TempFull(startIndex:startIndex+size(Temp,1)/4-1,:,:,:)=(Temp(1:4:size(Temp,1),presEnd:end,latUse,:)+...
    Temp(2:4:size(Temp,1),presEnd:end,latUse,:)+Temp(3:4:size(Temp,1),presEnd:end,latUse,:)...
    +Temp(4:4:size(Temp,1),presEnd:end,latUse,:))/4;
Geop=(nc_varget([zDir 'ERA5_z' fileAppend],'z'));
GeopFull(startIndex:startIndex+size(Temp,1)/4-1,:,:,:)=(Geop(1:4:size(Geop,1),presEnd:end,latUse,:)+...
    Geop(2:4:size(Geop,1),presEnd:end,latUse,:)+Geop(3:4:size(Geop,1),presEnd:end,latUse,:)...
    +Geop(4:4:size(Geop,1),presEnd:end,latUse,:))/4;


timeFull(startIndex:startIndex+size(Temp,1)/4-1)=time(1:4:length(time));
startIndex=startIndex+size(Temp,1)/4;

clear Geop Temp

end

lsmUse=reshape(repmat(land(:,latUse,:),length(timeFull)*length(pres),1,1),length(timeFull),length(pres),length(latUse),length(lon(2:end-1)));



lsmUse=cat(4,lsmUse(:,:,:,end),lsmUse,lsmUse(:,:,:,1));
Temp=cat(4,TempFull(:,:,:,end),TempFull,TempFull(:,:,:,1));%.*lsmUse;
Geop=cat(4,GeopFull(:,:,:,end),GeopFull,GeopFull(:,:,:,1));%.*lsmUse;

clear TempFull GeopFull

[times,P,Y,X]=ndgrid(timeFull,pres,lat,lon);
dx=X(:,:,:,3:end)-X(:,:,:,1:end-2);
dy=Y(:,:,3:end,:)-Y(:,:,1:end-2,:);
dp=diff(P,1,2);
dxkm=deg2rad(mod(dx,360)).*earthRad.*cosd(Y(:,:,:,2:end-1));
dykm=deg2rad(mod(dy,360))*earthRad;


S=(Temp*Cp+Geop);
S_bar=squeeze(mean(S,1));
dS_bar_dp=diff(S_bar,1,1)./squeeze(mean(dp,1));
S_tot=S_tot+S_bar;


eof(:,1) = lds(:,1)./sqrt(weights);
eof(:,2) = lds(:,2)./sqrt(weights);
eof(:,3) = lds(:,3)./sqrt(weights);
eof(:,4) = lds(:,4)./sqrt(weights);

Q1_O1 = repmat(eof(:,1),1,size(dS_bar_dp,2),size(dS_bar_dp,3)).*dS_bar_dp;
Q1_O2 = repmat(eof(:,2),1,size(dS_bar_dp,2),size(dS_bar_dp,3)).*dS_bar_dp; 

  
save([qDir 'ERA5_q1_year_mean' fileAppend(1:5) '.mat'],'Q1_O1','Q1_O2','lon','lat','presplot');



disp(years(yearsIndex));
%end
end
%end

S_bar_tot=S_tot./yearsIndex;
S_bar_tropics=squeeze(mean(mean(S_bar_tot,2),3));

dS_bar_tot_dp=diff(S_bar_tot,1,1)./squeeze(mean(dp,1));
dS_bar_tropics_dp=diff(S_bar_tropics,1,1)./squeeze(mean(mean(mean(dp,1),3),4));


Q1_O1_tot = repmat(eof(:,1),1,size(dS_bar_tot_dp,2),size(dS_bar_tot_dp,3)).*dS_bar_tot_dp;
Q1_O2_tot = repmat(eof(:,2),1,size(dS_bar_tot_dp,2),size(dS_bar_tot_dp,3)).*dS_bar_tot_dp; 
Q1_O3_tot = repmat(eof(:,3),1,size(dS_bar_tot_dp,2),size(dS_bar_tot_dp,3)).*dS_bar_tot_dp;
Q1_O4_tot = repmat(eof(:,4),1,size(dS_bar_tot_dp,2),size(dS_bar_tot_dp,3)).*dS_bar_tot_dp; 

Q1_O1_tropics = eof(:,1).*dS_bar_tropics_dp;
Q1_O2_tropics = eof(:,2).*dS_bar_tropics_dp; 
Q1_O3_tropics = eof(:,3).*dS_bar_tropics_dp;
Q1_O4_tropics = eof(:,4).*dS_bar_tropics_dp; 

save([qDir 'ERA5_q1_tot_mean.mat'],'Q1_O1_tot','Q1_O2_tot','Q1_O3_tot','Q1_O4_tot','lon','lat','presplot');
save([qDir 'ERA5_q1_tropics_mean.mat'],'Q1_O1_tropics','Q1_O2_tropics','Q1_O3_tropics','Q1_O4_tropics','presplot');
save(['ERA5_dsdp_tropics_mean.mat'],'lon','lat','preplot','dS_bar_tot_dp');
