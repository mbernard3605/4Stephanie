%function ERA5_PC() calculates the PCs for vertical motion by projecting the oceanic EOFs onto the vertical motion profiles
  %if running as a compiled program, this is needs to be turned off and these folders need to be added at compliation
if ~isdeployed
addpath('~/iris-home/matlab_scripts/mexcdf/mexnc');
addpath('~/iris-home/matlab_scripts/mexcdf/snctools');
end
%load('NCEP_EOFS.mat','omegaMean');
load('ERA5_dsdp_tropics_mean.mat');

%define constants
g=9.81;
earthRad=6.371e6;
Cp=1005;

%directories where all the variables and output are stored
baseDir='../Raw/ERA5/';
wDir = [baseDir 'w/'];
plotDir='../Plots/ERA5/';
qDir = [baseDir 'q1/'];


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
Q1=[];
%start loop through climatology
for yearsIndex = 1:length(years)
    startIndex=1;
    outFile=[qDir 'Q1_wtg_' num2str(years(yearsIndex)) '.mat'];
%create netcdf outfile if one does not exist
if ~exist(outFile)
    if ~mod(years(yearsIndex),4)
        baseFile='../baseleap_ERA5.nc';
        days=366;
        OmegaFull=zeros(days,27,19,144);
        timeFull=zeros(days,1);
    else
        baseFile='../base_ERA5.nc';
        days=365;
        OmegaFull=zeros(days,27,19,144);
        timeFull=zeros(days,1);
    end
copyfile(baseFile,outFile);
end
%loop through each month this year
for monthIndex = 1:length(months);

            fileAppend=sprintf('_%d_%d.nc',years(yearsIndex),months(monthIndex));
time = nc_varget([wDir 'ERA5_w' fileAppend],'time');
%if this is the first go around we define constants
if ~started

lon=nc_varget([wDir 'ERA5_w' fileAppend],'longitude');
lat=nc_varget([wDir 'ERA5_w' fileAppend],'latitude');
latUse = find(lat<=22.5 & lat>=-22.5); 
lon = [lon(end); lon; lon(1)];
lat=lat(latUse);
level=nc_varget([wDir 'ERA5_w' fileAppend],'level')*100;
presEnd=find(level==10000);
dp=diff(level);
weights=dp./sum(dp);

pres=level(presEnd:end);
started=1;


[xplot,yplot]=meshgrid(lon(2:end-1),lat(2:end-1));
presplot=(pres(2:end)+pres(1:end-1))/2;

end



%load the vertical velocity
Omega=nc_varget([wDir 'ERA5_w' fileAppend],'w');
OmegaFull(startIndex:startIndex+size(Omega,1)/4-1,:,:,:)=(Omega(1:4:size(Omega,1),presEnd:end,latUse,:)+...
    Omega(2:4:size(Omega,1),presEnd:end,latUse,:)+Omega(3:4:size(Omega,1),presEnd:end,latUse,:)...
    +Omega(4:4:size(Omega,1),presEnd:end,latUse,:))/4;


timeFull(startIndex:startIndex+size(Omega,1)/4-1)=time(1:4:length(time));
startIndex=startIndex+size(Omega,1)/4;

clear Omega

end
%creat land sea mask if it is desired
lsmUse=reshape(repmat(land(:,latUse,:),length(timeFull)*length(pres),1,1),length(timeFull),length(pres),length(latUse),length(lon(2:end-1)));



lsmUse=cat(4,lsmUse(:,:,:,end),lsmUse,lsmUse(:,:,:,1));
Omega=cat(4,OmegaFull(:,:,:,end),OmegaFull,OmegaFull(:,:,:,1));%.*lsmUse;

clear OmegaFull

%decompose the vertical motion into the first two EOFs

%save the years's PCs to the matrix of all the PCs

%save this years PCs
%save([baseDir 'PC/' 'ERA5_PC' fileAppend(1:5) '.mat'],'O1','O2');
Q1=Omega.*permute(repmat(dS_bar_tot_dp,[1,1,1,size(Omega,1)]),[4 1 2 3]);

save(outFile,'Q1');
disp(years(yearsIndex));

%end
end
%end

%save the large matrices of the climatology of PCs



