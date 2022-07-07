clear;

%% Post-process Multi-surface Run

% Specify Run to Process
runDate='20211115_ADJ_M3/'; %choose run to process here
runname= 'multi-adj-ekh-alb-2007'; % change run name too

% ICEMELT Output Directory
outDirectory='output/';
path2output=[outDirectory runDate];

% Directories to process
surfrun=       [path2output 'XXX_' runname, '_0.0000500_0.0020_09/' ];
basinwallrun=  [path2output 'XXX_' runname '-bwall_0.0010000_0.0020_09/'];
basinfloorrun= [path2output 'XXX_' runname '-bfloor_0.0010000_0.0020_09/'];
cliffrun=      [path2output 'XCL_' runname '-cliff_0.0001000_0.0020_09/'];

% surfrun=       [path2output runname, '/' ];
% basinwallrun=  [path2output runname '-bwall/'];
% basinfloorrun= [path2output runname '-bfloor/'];
% cliffrun=      [path2output runname '-cliff/'];

%% Script details - what to process, where some data files are located.

% Meltwater input data
modelinputdir  = 'input/';
streamdatafile = 'input/TV_stream_data.mat';
dataqualfile = 'input/TV_stream_data_qual.mat';

docliff=1;  % option to include analysis of cliff model or not.
dobasins=1;  % option to include analysis of basins or not.

adjustBasins = 1;

years2run = (1995:2012);

yearIndices = years2run-1994;

ndays= datenum([years2run(end)+1 1 31]) - datenum([1995 6 30]);
% Model runs to 02/01/2012

nx=200; ny=140;
dx=250; dy=dx;

basins2process=(1:99); % whole valley
%basins2process=(41);  % CAA

%%
%[basins,Ref] = ImportAsciiRaster([ modelinputdir 'tv_basins_surface.txt']);
[basins,Ref] = ImportAsciiRaster([ modelinputdir 'tv_basins_surface_jmc_ekh.txt']);

% Note: Priscu stream's contributing area needs to be adjusted
% Priscu stream fed by multiple glaciers (use code Basin 21)
i=find( (basins==22) | (basins==23) );    basins(i)=21;

% Note: Include some of Commonwealth (73) into Aiken 
% streamflow (e.g. outflow from Many Glacier Pond)
i=find( (basins==62) | (basins==73) );    basins(i)=62;

% Stream Basin Key
% streamkey=[62 41 45 74 66 50 43 65 34 61 29 71 21 10 63];
streamkey=[62 41 45 74 66 50 43 65 34 61 29 71 21 10 64];

streamname={'aiken';'andersen';'canada';'common';'crescent';...
    'delta';'green';'harnish';'house';'huey';...
    'lawson';'lostseal';'priscu';'santafe';'vonguerard' };


for thebasin=1:15
    numcells(thebasin)=sum(sum(basins==streamkey(thebasin)));
end

% Read Stream Data
load(streamdatafile);
load(dataqualfile);
streams.dates=streams.dates(1:ndays);  % starts at 1 jul 1995
streams.TDQ=streams.TDQ(1:ndays,:);
streams.TDQ=streams.TDQ * 0.001;  %L to m^3 conversion
streams.ErrUpr=streams.TDQ+(streams.TDQ.*error.Error(1:ndays,:));
streams.ErrLwr=streams.TDQ-(streams.TDQ.*error.Error(1:ndays,:));

basin_inds=[];
for i=1:length(basins2process)
    basin_inds=[basin_inds; find(streamkey==basins2process(i))]; %#ok<*AGROW>
end

% Read y coords for this run
M=dlmread([surfrun 'ycrds.out']);
subdepths=M(2:length(M)-1,1);
%THIS ISNT RIGHT BUT SHOULD BE CLOSE
deltaz(1)=subdepths(1)-0;
for i=2:length(subdepths)
    deltaz(i)=subdepths(i)-subdepths(i-1);
end

density=0.87;


%% Figure out areas of each subregion

fprintf('Calculating area of each domain... \n');

% Load the ASCII grids file that says where the basins are.
[basin_ridge_area, Ref] = ImportAsciiRaster([ modelinputdir 'tv_bridge_number_of_2mgridcells_per_250mgridcell.txt']); %#ok<*ASGLU>
[basin_walls_area, Ref] = ImportAsciiRaster([ modelinputdir 'tv_bwalls_number_of_2mgridcells_per_250mgridcell.txt']);
[basin_floor_area, Ref] = ImportAsciiRaster([ modelinputdir 'tv_bfloor_number_of_2mgridcells_per_250mgridcell.txt']);
%Turn NaN into 0.
%cleaner=find(isnan(basin_all_area));   basin_all_area(cleaner)=0;  % change NaN to 0
cleaner=find(isnan(basin_ridge_area));   basin_ridge_area(cleaner)=0;  %#ok<*FNDSB> % change NaN to 0
cleaner=find(isnan(basin_walls_area));   basin_walls_area(cleaner)=0;  % change NaN to 0
cleaner=find(isnan(basin_floor_area));   basin_floor_area(cleaner)=0;  % change NaN to 0

% Convert from number of 2m grid cells to an area
if docliff==1
    [basinsCliff, Ref] = ImportAsciiRaster([ modelinputdir '/tv_basins_cliff.txt']);
    [cliffLength, Ref] = ImportAsciiRaster([ modelinputdir 'cliff_length.txt']);
    cleaner=find(cliffLength(:)==0);
    cliffLength(cleaner)=NaN;
    cliffHeight=25; %maybe need to refine this!
end

basin_wall_area_factor = 1.2;  % scale to adjust area of basin walls - (increase basin walls by 20% due to high slopes)

%basin_all_area=basin_all_area*2*2;
basin_ridge_area=basin_ridge_area*2*2;
basin_walls_area=basin_walls_area*2*2 * basin_wall_area_factor;  % adjust basin wall area by a factor
basin_floor_area=basin_floor_area*2*2;

% Make the uniform surface a combo of uniform surfaces and basin ridges
uniformsurface_area=ones(size(basin_ridge_area)) * 250 * 250;
% start off my assuming a full size grid cell everywhere.
uniformsurface_area=uniformsurface_area - basin_walls_area/basin_wall_area_factor - basin_floor_area;
% remove the basin areas, but keep the basin ridges

% MJH: NOTE THE UNIFORM SURFACE USE THE WHOLE GRID CELL
% EVEN IF THE GLACIER OUTLINE SPLITS THE CELL!
% THIS APPEARS TO OVERESTIMATE MOST BASINS

% cryo hole area
cryohole_area=uniformsurface_area*0.043;  %4.3% is for Canada
uniformsurface_area = uniformsurface_area - cryohole_area;  % dont double count these areas

if docliff==1
    cliff_area = cliffLength .* cliffHeight;
end

noncliff_area = uniformsurface_area + basin_ridge_area + basin_walls_area + basin_floor_area + cryohole_area;
%% --Area Information--

% b = 34;  % pick the basin to analyze
%
% i=find(basins==b);
% ic=find(basinsCliff==b);
%
% a = sum(noncliff_area(i)) + sum(cliff_area(ic))
%
% sum(basin_ridge_area(i)) / a
% sum(basin_walls_area(i)) / a
% sum(basin_floor_area(i)) / a
% sum(cryohole_area(i)) / a
% sum(cliff_area(ic)) / a

%% Create placeholders for the new data

nyears = years2run(end) - years2run(1) + 1;

% Smooth and Cryoholes
uniformsurface.surfTDQ=zeros(size(streams.TDQ));
uniformsurface.subTDQ=zeros(size(streams.TDQ));
cellmeltsummer=zeros(nyears,nx,ny)*NaN;
cellablsummer=zeros(nyears,nx,ny)*NaN;
cellsubmeltsummer=zeros(nyears,nx,ny)*NaN;
cellsubmeltSurfacedsummer=zeros(nyears,nx,ny)*NaN;

% Basin Walls
basinwall.surfTDQ=zeros(size(streams.TDQ));
basinwall.subTDQ=zeros(size(streams.TDQ));
basinwall.cellmeltsummer=zeros(nyears,nx,ny)*NaN;
basinwall.cellablsummer=zeros(nyears,nx,ny)*NaN;
basinwall.cellsubmeltsummer=zeros(nyears,nx,ny)*NaN;
basinwall.cellsubmeltSurfacedsummer=zeros(nyears,nx,ny)*NaN;

%Basin Floors
basinfloor.surfTDQ=zeros(size(streams.TDQ));
basinfloor.subTDQ=zeros(size(streams.TDQ));
basinfloor.cellmeltsummer=zeros(nyears,nx,ny)*NaN;
basinfloor.cellablsummer=zeros(nyears,nx,ny)*NaN;
basinfloor.cellsubmeltsummer=zeros(nyears,nx,ny)*NaN;
basinfloor.cellsubmeltSurfacedsummer=zeros(nyears,nx,ny)*NaN;

% Cliffs
cliffTDQ=zeros(size(streams.TDQ));
cliffs.surfTDQ=zeros(size(streams.TDQ));
cliffs.subTDQ=zeros(size(streams.TDQ));
cliffarea=zeros(size(streamkey));
cliffcells=zeros(size(streamkey));
cliffMnAnnlMelt=zeros(size(streamkey));
totalcliffarea=zeros(1,15);
cliff.cellmeltsummer=zeros(length(yearIndices),nx,ny)*NaN;
cliff.cellsubmeltsummer=zeros(length(yearIndices),nx,ny)*NaN;

% Totals
cryoholeTDQ=zeros(size(streams.TDQ));
surfTDQ=zeros(size(streams.TDQ));
subTDQ=zeros(size(streams.TDQ));
modelTDQ=zeros(size(streams.TDQ)); %#ok<*PREALL>

% Other
surfarea=zeros(size(streamkey));
surfcells=zeros(size(streamkey));
surfMnAnnlMelt=zeros(size(streamkey));

%% Read Melt Files            SMOOTH

fprintf('Running smooth surfaces domain... \n');

for i=1:nx
    for j=1:ny
        
        % Surface model output
        if ismember(basins(ny-j+1,i),basins2process)  %y index is backward from the ascii file
            iname=num2str(i,'%3.3d');
            jname=num2str(j,'%3.3d');
            fname=[surfrun ,iname,jname,'.ablation.daily'];
            fid=fopen(fname,'r');
            M=fread(fid,[7,Inf],'float32')'; %daily melt, ablation
            fclose(fid);
            
            % Calculate subsurface melt contribution to surface ablation
            % read in density loss for whole season
            endofsummerdensity=dlmread([surfrun '/',iname,jname,'.densityprofile']);
            season=2;
            
            for yy=years2run
                % Calculate seasonal totals for all cells every year
                seasonstart=datenum([yy 11 1])-datenum([1995 6 30]);
                if yy==years2run(length(years2run))
                    seasonend=datenum([yy+1 1 31])-datenum([1995 6 30]);
                else
                    seasonend=datenum([yy+1 2 28])-datenum([1995 6 30]);
                end
                cellmeltsummer(yy-1994,i,j)=sum(M(seasonstart:seasonend,1));
                cellablsummer(yy-1994,i,j)=sum(M(seasonstart:seasonend,2));
                cellsubmeltsummer(yy-1994,i,j)=sum(M(seasonstart:seasonend,3));
                
                % Find depth to query the subsurface densities over
                affectedlayers=interp1(subdepths,(1:length(subdepths)),cellablsummer(yy-1994,i,j)/100);
                temp=0;
                if (affectedlayers>29)
                    affectedlayers=29;
                end
                for ii=1:floor(affectedlayers)
                    temp=temp+deltaz(ii)*(density*1000-endofsummerdensity(yy-1994,ii))/1000;
                end
                ii=ii+1;
                temp=temp+deltaz(ii)*(density*1000-endofsummerdensity(yy-1994,ii))/1000 * (affectedlayers-floor(affectedlayers));
                % Fractional part of final layer
                cellsubmeltSurfacedsummer(yy-1994,i,j)=temp*100;  % save it
            end
            
            %Apportion this melt to the proper basin
            basinnum=find(streamkey==basins(ny-j+1,i));
            if length(basinnum==1)  % Not all cells belong to basins
                %Smooth Surface Melt
                uniformsurface.surfTDQ(:,basinnum)=uniformsurface.surfTDQ(1:ndays,basinnum) + uniformsurface_area(ny-j+1,i)*M(1:ndays,1)/100;
                % Smooth Sub Melt
                uniformsurface.subTDQ(:,basinnum)=uniformsurface.subTDQ(1:ndays,basinnum) + uniformsurface_area(ny-j+1,i)*M(1:ndays,3)/100;
                % Cryoholes
                cryoholeTDQ(:,basinnum)=cryoholeTDQ(1:ndays,basinnum) + cryohole_area(ny-j+1,i) * M(1:ndays,2)/100;
                for yy=years2run
                    seasonstart=datenum([yy 11 1])-datenum([1995 6 30]);
                    if yy==years2run(length(years2run))
                        seasonend=datenum([yy+1 1 31])-datenum([1995 6 30]);
                    else
                        seasonend=datenum([yy+1 2 28])-datenum([1995 6 30]);
                    end
                    cryoholeTDQ(seasonstart:seasonend,basinnum)=cryoholeTDQ(seasonstart:seasonend,basinnum) + cryohole_area(ny-j+1,i) * cellsubmeltSurfacedsummer(yy-1994,i,j)/100/(seasonend-seasonstart+1); % just spread it evenly over the season...
                end
            end
        end
    end
end

cryoholeTDQ=max(cryoholeTDQ,cryoholeTDQ*0);  % make sure this is not negative
surfMnAnnlMelt=surfMnAnnlMelt./surfcells;
cellcryoholesummer = max(cellablsummer, cellablsummer*0.0);

%% Read Melt Files            BASIN WALL

fprintf('Running basin wall domain... \n');

if dobasins == 1
    
    for i=1:nx
        for j=1:ny
            
            if adjustBasins == 1
                if basins(ny-j+1,i) == 10 || basins(ny-j+1,i) == 15 || basins(ny-j+1,i) == 19
                    adjFactor = 1.25;
                elseif basins(ny-j+1,i) == 29 || basins(ny-j+1,i) == 21
                    adjFactor = 3;
                elseif basins(ny-j+1,i) > 40 && basins(ny-j+1,i) < 46
                    adjFactor = 0.4;
                else
                    adjFactor = 1;
                end
            else
                adjFactor = 1;
            end
            
            % Surface model output
            if ismember(basins(ny-j+1,i),basins2process)  %y index is backward from the ascii file
                iname=num2str(i,'%3.3d');
                jname=num2str(j,'%3.3d');
                fname=[basinwallrun ,iname,jname,'.ablation.daily'];
                fid=fopen(fname,'r');
                M=fread(fid,[7,Inf],'float32')'; %daily melt, ablation
                fclose(fid);
                
                % Calculate subsurface melt contribution to surface ablation
                % read in density loss for whole season
                endofsummerdensity=dlmread([surfrun '/',iname,jname,'.densityprofile']);
                season=2;
                
                for yy=years2run
                    % Calculate seasonal totals for each cell
                    seasonstart=datenum([yy 11 1])-datenum([1995 6 30]);
                    if yy==years2run(length(years2run))
                        seasonend=datenum([yy+1 1 31])-datenum([1995 6 30]);
                    else
                        seasonend=datenum([yy+1 2 28])-datenum([1995 6 30]);
                    end
                    basinwall.cellmeltsummer(yy-1994,i,j)=sum(M(seasonstart:seasonend,1));
                    basinwall.cellablsummer(yy-1994,i,j)=sum(M(seasonstart:seasonend,2));
                    basinwall.cellsubmeltsummer(yy-1994,i,j)=sum(M(seasonstart:seasonend,3));
                    
                    % Find depth to query the subsurface densities over
                    affectedlayers=interp1(subdepths,[1:length(subdepths)],cellablsummer(yy-1994,i,j)/100);
                    temp=0;
                    if (affectedlayers>29); affectedlayers=29; end
                    for ii=1:floor(affectedlayers)
                        temp=temp+deltaz(ii)*(density*1000-endofsummerdensity(yy-1994,ii))/1000;
                    end
                    ii=ii+1;
                    temp=temp+deltaz(ii)*(density*1000-endofsummerdensity(yy-1994,ii))/1000 * (affectedlayers-floor(affectedlayers));  %fractional part of final layer
                    basinwall.cellsubmeltSurfacedsummer(yy-1994,i,j)=temp*100;  % save it
                end
                
                % Apportion this melt to the proper basin
                basinnum=find(streamkey==basins(ny-j+1,i));
                if length(basinnum==1)  %&&(basin_regions(ny-j+1,i)>=0)  %Not all cells belong to basins, check to make sure this is a basin-y part of the drainage basin
                    basinwall.surfTDQ(1:ndays,basinnum)=basinwall.surfTDQ(1:ndays,basinnum)+ adjFactor*basin_walls_area(ny-j+1,i)*M(1:ndays,1)/100; % Correct for grid size and units - i want cubic meters
                    basinwall.subTDQ(1:ndays,basinnum)=basinwall.subTDQ(1:ndays,basinnum)+ adjFactor*basin_walls_area(ny-j+1,i)*M(1:ndays,3)/100; %subsurf melt is in m :(
                end
            end
        end
    end
    
end

%% Read Melt Files            BASIN FLOOR

fprintf('Running basin floor domain... \n');

if dobasins == 1
    
    for i=1:nx
        for j=1:ny
            
            if adjustBasins == 1
                if basins(ny-j+1,i) == 10 || basins(ny-j+1,i) == 15 || basins(ny-j+1,i) == 19
                    adjFactor = 1.25;
                elseif basins(ny-j+1,i) == 29 || basins(ny-j+1,i) == 21
                    adjFactor = 3;
                elseif basins(ny-j+1,i) > 40 && basins(ny-j+1,i) < 46
                    adjFactor = 0.4;
                else
                    adjFactor = 1;
                end
            else
                adjFactor = 1;
            end
            
            % Surface model output
            if ismember(basins(ny-j+1,i),basins2process)  %y index is backward from the ascii file
                iname=num2str(i,'%3.3d');
                jname=num2str(j,'%3.3d');
                fname=[basinfloorrun ,iname,jname,'.ablation.daily'];
                fid=fopen(fname,'r');
                M=fread(fid,[7,Inf],'float32')'; %daily melt, ablation
                fclose(fid);
                
                % Calculate subsurface melt contribution to surface ablation
                % read in density loss for whole season
                endofsummerdensity=dlmread([surfrun '/',iname,jname,'.densityprofile']); % Note JMC: Turned this back on?!?
                season=2;
                for yy=years2run
                    % Calculate seasonal totals for each cell.
                    seasonstart=datenum([yy 11 1])-datenum([1995 6 30]);
                    if yy==years2run(length(years2run))
                        seasonend=datenum([yy+1 1 31])-datenum([1995 6 30]);
                    else
                        seasonend=datenum([yy+1 2 28])-datenum([1995 6 30]);
                    end
                    basinfloor.cellmeltsummer(yy-1994,i,j)=sum(M(seasonstart:seasonend,1));
                    basinfloor.cellablsummer(yy-1994,i,j)=sum(M(seasonstart:seasonend,2));
                    basinfloor.cellsubmeltsummer(yy-1994,i,j)=sum(M(seasonstart:seasonend,3));
                    
                    % Find depth to query the subsurface densities over
                    affectedlayers=interp1(subdepths,[1:length(subdepths)],cellablsummer(yy-1994,i,j)/100);
                    temp=0;
                    if (affectedlayers>29); affectedlayers=29; end
                    for ii=1:floor(affectedlayers)
                        temp=temp+deltaz(ii)*(density*1000-endofsummerdensity(yy-1994,ii))/1000;
                    end
                    ii=ii+1;
                    temp=temp+deltaz(ii)*(density*1000-endofsummerdensity(yy-1994,ii))/1000 * (affectedlayers-floor(affectedlayers));  %fractional part of final layer
                    basinfloor.cellsubmeltSurfacedsummer(yy-1994,i,j)=temp*100;  % save it
                end
                
                % Apportion this melt to the proper basin
                basinnum=find(streamkey==basins(ny-j+1,i));
                if length(basinnum==1)% Not all cells belong to basins, check to make sure this is a basin-y part of the drainage basin
                    basinfloor.surfTDQ(1:ndays,basinnum)=basinfloor.surfTDQ(1:ndays,basinnum)+ adjFactor*basin_floor_area(ny-j+1,i)*M(1:ndays,1)/100; % Correct for grid size and units - i want cubic meters
                    basinfloor.subTDQ(1:ndays,basinnum)=basinfloor.subTDQ(1:ndays,basinnum)+ adjFactor*basin_floor_area(ny-j+1,i)*M(1:ndays,3)/100; %subsurf melt is in m :(
                end
            end
        end
    end
    
end

%% Do cliff model            CLIFF

fprintf('Running cliff domain... \n');

if docliff==1
    
    for i=1:nx
        
        for j=1:ny
            
            % Cliff model output
            if ismember(basinsCliff(ny-j+1,i),basins2process)  %y index is backward from the ascii file
                iname=num2str(i,'%3.3d');
                jname=num2str(j,'%3.3d');
                fname=[cliffrun,iname,jname,'.ablation.daily'];
                fid=fopen(fname,'r');
                M=fread(fid,[7,Inf],'float32')'; %daily melt, ablation
                fclose(fid);
                
                season = 2;
                for yy=years2run
                    seasonstart=datenum([yy 11 1])-datenum([1995 6 30]);
                    if yy==years2run(length(years2run))
                        seasonend=datenum([yy+1 1 31])-datenum([1995 6 30]);
                    else
                        seasonend=datenum([yy+1 2 28])-datenum([1995 6 30]);
                    end
                    cliff.cellmeltsummer(yy-1994,i,j)=sum(M(seasonstart:seasonend,1));
                    cliff.cellablsummer(yy-1994,i,j)=sum(M(seasonstart:seasonend,2));
                    cliff.cellsubmeltsummer(yy-1994,i,j)=sum(M(seasonstart:seasonend,3));
                end %for year
                
                %Apportion this melt to the proper basin
                basinnum=find(streamkey==basinsCliff(ny-j+1,i));
                if length(basinnum==1)  %Not all cells belong to basins
                    cliffs.surfTDQ(1:ndays,basinnum)=cliffs.surfTDQ(1:ndays,basinnum)+M(1:ndays,1)/100*cliffLength(ny-j+1,i)*cliffHeight; % Correct for cliff height and units
                    cliffs.subTDQ(1:ndays,basinnum)=cliffs.subTDQ(1:ndays,basinnum)+M(1:ndays,3)/100*cliffLength(ny-j+1,i)*cliffHeight; %subsurf melt
                    cliffarea(basinnum)=cliffarea(basinnum) + cliffLength(ny-j+1,i)*cliffHeight;
                    cliffcells(basinnum)=cliffcells(basinnum)+1;
                    cliffMnAnnlMelt(basinnum)=cliffMnAnnlMelt(basinnum)+sum(M(:,1))/13;
                    totalcliffarea(basinnum)=totalcliffarea(basinnum)+cliffLength(ny-j+1,i)*cliffHeight;
                end
            end
        end
    end
    cliffMnAnnlMelt=cliffMnAnnlMelt./cliffcells;
    cliffTDQ=cliffs.surfTDQ+cliffs.subTDQ;  % keep this old variable , but add in both parts
end

%% Calculate Totals

fprintf('Calculating totals... \n');

% Calculate domain melt season totals
for y=yearIndices
    %Summer melt
    domainsummermelt(y,1)=nansum(nansum(nansum(  squeeze(cellmeltsummer(y,:,:)/100)'.* flipud(uniformsurface_area)   )));
    domainsummermelt(y,2)=nansum(nansum(nansum(  squeeze(cellsubmeltsummer(y,:,:)/100)'.* flipud(uniformsurface_area)   )));
    domainsummermelt(y,3)=nansum(nansum(nansum(  squeeze(basinwall.cellmeltsummer(y,:,:)/100)'.* flipud(basin_walls_area)   )));
    domainsummermelt(y,4)=nansum(nansum(nansum(  squeeze(basinwall.cellsubmeltsummer(y,:,:)/100)'.* flipud(basin_walls_area)   )));
    domainsummermelt(y,5)=nansum(nansum(nansum(  squeeze(basinfloor.cellmeltsummer(y,:,:)/100)'.* flipud(basin_floor_area)   )));
    domainsummermelt(y,6)=nansum(nansum(nansum(  squeeze(basinfloor.cellsubmeltsummer(y,:,:)/100)'.* flipud(basin_floor_area)   )));
    %domainsummermelt(y,7)=sum(sum(sum(  cliff.cellmeltsummer(y,:,:) )));  % in cubic m
    %domainsummermelt(y,8)=sum(sum(sum(  cliff.cellsubmeltsummer(y,:,:) ))); % in cubic m
end

%Calculate Total Volume for each season
modelTDQsubonly=uniformsurface.subTDQ + basinfloor.subTDQ + basinwall.subTDQ; %JMC Changed these!!
modelTDQsurfonly=uniformsurface.surfTDQ + basinfloor.surfTDQ + basinwall.surfTDQ;
modelTDQ = modelTDQsubonly + modelTDQsurfonly + cryoholeTDQ + cliffTDQ;

%modelTDQ=uniformsurface.surfTDQ+cliffTDQ+basinfloor.surfTDQ+basinwall.surfTDQ + uniformsurface.subTDQ + basinfloor.subTDQ+basinwall.subTDQ + cryoholeTDQ;

% Loop through each basin
for b=[1:length(streamkey)]
    for yr=years2run
        
        % Set season start and end
        strt=find(streams.dates==datenum([yr 7 1]));
        fnsh=find(streams.dates==datenum([yr+1 6 30]));
        
        % Set summer start and end
        sumrstart=find(streams.dates==datenum([yr 10 1]));
        sumrfnsh=find(streams.dates==datenum([yr+1 3 31]));
        
        % Use last date available
        if (yr==years2run(end))
            fnsh=length(modelTDQ);
            sumrfnsh=length(modelTDQ);
        end
        
        seasonDischarge = nansum(streams.TDQ(strt:fnsh,b));
        
        if size(seasonDischarge) == 1
            streamYrVol(yr-1994,b)=seasonDischarge;
        else
            streamYrVol(yr-1994,b)=NaN;  % SHOULD THIS BE NAN?
        end
        measdays=find( ~isnan(streams.TDQ(:,b)) & ([1:length(streams.dates)]'>=strt) & ([1:length(streams.dates)]'<=fnsh) );  % note days that have data only
        
        % Model totals
        modelYrVol(yr-1994,b)=sum(modelTDQ(strt:fnsh,b));
        modelSmVol(yr-1994,b)=sum(modelTDQ(sumrstart:sumrfnsh,b));
        
        % Surface and Subsurface totals
        surfYrVol(yr-1994,b)=sum(modelTDQsurfonly(strt:fnsh,b));
        subYrVol(yr-1994,b)=sum(subTDQ(strt:fnsh,b));
        
        % Surface Type Totals
        cliffYrVol(yr-1994,b)=sum(cliffTDQ(strt:fnsh,b));
        uniformsurfaceYrVol(yr-1994,b)=sum(uniformsurface.surfTDQ(strt:fnsh,b) + uniformsurface.subTDQ(strt:fnsh,b));
        basinsYrVol(yr-1994,b)=sum(basinfloor.surfTDQ(strt:fnsh,b) + basinfloor.subTDQ(strt:fnsh,b) + basinwall.surfTDQ(strt:fnsh,b) + basinwall.subTDQ(strt:fnsh,b));
        cryoholeYrVol(yr-1994,b)=sum(cryoholeTDQ(strt:fnsh,b));
        
        % Only use days that have obs
        modelYrVolMeas(yr-1994,b)=sum(modelTDQ(measdays,b));
        % surfYrVol(yr-1994,b)=sum(modelTDQsurfonly(measdays,b));
    end
end

%% SAVE OUTPUT

fprintf('Saving output...\n');

if adjustBasins == 1
    %runDate='20200225_BASINS_M4/'; %choose run to process here
    runname = ['basin-' runname];
end

% Data Output Directory
outDirectory='processed/';

fullPathData = fullfile(outDirectory, runDate);

if ~exist(fullPathData, 'dir')
    
    % Folder does not exist so create it.
    mkdir (fullPathData);
    
    fileName = [runname '-streams.mat'];
    path2save=[outDirectory runDate fileName];
    
    % Save output data
    save(path2save, 'streamkey', 'streamname', 'streams', 'streamYrVol', 'modelTDQ', 'modelYrVol', 'modelYrVolMeas', ...
    'cliffYrVol', 'cryoholeYrVol', 'uniformsurfaceYrVol', 'basinsYrVol', 'surfYrVol', 'modelTDQsurfonly', ...
    'cliffTDQ', 'cryoholeTDQ', 'uniformsurface', 'basinfloor', 'basinwall', 'modelSmVol');
else

    fileName = [runname '-streams.mat'];
    path2save=[outDirectory runDate fileName];
    
    % Save output data
    save(path2save, 'streamkey', 'streamname', 'streams', 'streamYrVol', 'modelTDQ', 'modelYrVol', 'modelYrVolMeas', ...
        'cliffYrVol', 'cryoholeYrVol', 'uniformsurfaceYrVol', 'basinsYrVol', 'surfYrVol', 'modelTDQsurfonly', ...
        'cliffTDQ', 'cryoholeTDQ', 'uniformsurface', 'basinfloor', 'basinwall', 'modelSmVol');
end
%% End