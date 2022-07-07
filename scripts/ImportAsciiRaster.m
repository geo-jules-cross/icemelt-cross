function [Z R] = ImportAsciiRaster(varargin)
% [Z R] = ImportAsciiRaster(...)
% 
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%     %                                                                 %
%     %                 Produced by Giuliano Langella                   %
%     %                   e-mail:gyuliano@libero.it                     %
%     %                           March 2008                            %
%     %                                                                 %
%     %                 Last Updated: 20 April, 2010                    %
%     %                                                                 %
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
%
%
%   -----   TODO   -----
%   1.\ The function should take care of xll-center and yll-center when
%       they are found instead of -corner conditions.
% 
%   -----  SYNTAX  -----
%   Z = ImportAsciiRaster;                  ==> nodata:NaN; ref:'h'; type:'s'; FilePath:uigetfile
%   [Z h] = ImportAsciiRaster;              ==> nodata:NaN; ref:'h'; type:'s'; FilePath:uigetfile
%   [Z R] = ImportAsciiRaster(varargin);    ==> see OPTIONS
%
%
%   -----  OPTIONS  -----
%  +'nodata'(optional) [numeric] = output No-Data Value.
%   Default output NoData is NaN. When specified it substitutes the input
%   no-data value, which is stored in 'NODATA_value' from header.
%
%  +'ref'(optional) [string]: Geospatial Reference Matrix.
%   --> 'h' (default) for header matrix (as the .hdr file in Arc-Info Binary Raster).
%   --> 'r' for spatial-referencing matrix R (Mapping Toolbox);
%
%  +'type'(optional) [string]: class type of Z output in multiselection case only.
%   --> 's' (default) gets a structure array with field name equal to
%       imported filename without extension;
%   --> 'd' gets a double array with size (nrows, ncols, N. files); it is
%       useful in case of several layers with same reference matrix.
%       NOTE: be careful to multi-select grids with same spatial extent.
%
%  +'FilePath'(optional) [char]: Complete path + filename + file extension.
%       If it is not given the uigetfile will ask for file(s).
%
%   -----  DESCRIPTION  -----
% This function imports Arc-Info ASCII Raster (*.asc, *.txt) files. You are
% allowed to import single or multiple files by multi-selection (CRTL +
% Mouse_Left_Click). In the latter case by default a structure array is
% created in order to store the several rasters; items of the structure
% array are called as the imported raster. A facultative option is now
% added to create a stack of double 3-D array instead of a structure array
% in case more grid files are selected. This might be useful to handle for
% instance multi-temporal remote sensing datasets (such as NDVI), or
% multiple realizations of Geostatistical Simulation. You can create a
% stack with dimension compatibly with you computer memory.
% 
%   -----  EXAMPLE_1  -----
%   [Z R] = ImportAsciiRaster('r');           % import DEM and Aspect
%   Exaggeration_Factor = 1.6;                % use an exaggeration factor
%   figure(gcf); mapshow(Z.DEM*Exaggeration_Factor,R.DEM,'DisplayType','surface');
%   colormap(demcmap(Z.DEM))
%   view(3); grid off; axis off
%   axis on; xlabel('Easting'); ylabel('Northing'); zlabel('Elevation')
%
%   -----  EXAMPLE_2  -----
%   [Z R] = ImportAsciiRaster(NaN, 'r', 'd');       %import 15 Geostatistical Simulations
%   M = mean(Z,3);      %compute mean along 3rd dimension
%   figure(gcf); h=pcolor(M); axis off;             %pcolor plot
%   set(h,'LineStyle','none')
%   colormap(jet), colorbar('Location', 'SouthOutside')
%   saveas(gcf, 'clay_simulation.eps', 'psc2')      %save to color EPS for LaTex
%

%% MAIN

% Read varargins [given by user]
FilePath = cell(0);
FileName = cell(0);
for i = 1:size(varargin,2)
    if iscell(varargin{i})
        FilePath = varargin{i};
        for i = 1:length(FilePath)
            FileName{end+1} = path2file(FilePath{i});
        end
    else
        switch varargin{i}
            case{'h','r'}                                       % header
                ref = varargin{i};
            case{'s','d'}                                       % var type
                type = varargin{i};
            otherwise        
                if isnan(varargin{i}) | isnumeric(varargin{i})  % nodata
                    nodata = varargin{i};
                elseif ischar(varargin{i})                      % filepath
                    FilePath{end+1} = varargin{i};
                    FileName{end+1} = path2file(FilePath{end});
                    if ~exist(FilePath{end},'file')
                        error(['File ',varargin{i},' does not exist!'])
                    end
                end
        end
    end
end

% get file(s) [use UI if any file is given as varargin]
if ~sum(size(FilePath))
    % UIGETFILE - .asc/.txt - Multiselect on
    [FileName, PathName] = uigetfile({'*.asc','Arc-Info ASCII Raster (*.asc)'; ...
       '*.txt','Text ASCII Raster (*.txt)'; ...
       '*.*',  'All Files (*.*)'}, ...
       'Pick Raster(s)', ...
       pwd, ...
       'MultiSelect', 'on');
    % FILEPATH
    if iscell(FileName)
        for NumFile = 1:length(FileName)
            FilePath{NumFile} = strcat(PathName, FileName(NumFile));
        end
    else
        FilePath{1} = strcat(PathName, FileName);
        FileName = {FileName};
    end
end

% set undefined varargins
if ~exist('nodata', 'var'),     nodata  = NaN;      end     % default
if ~exist('ref', 'var'),        ref     = 'h';      end     % default
if ~exist('type', 'var'),       type    = 's';      end     % default
if length(FileName)==1,         type    =  '';      end     % double

%'''''''''''''''''
% clc
%fprintf('\nReading...\n')
%'''''''''''''''''

% START LOOP FOR EACH FILE TO BE IMPORTED
for NumFile = 1:size(FilePath,2)

    %fprintf('\t-->[%s]\n', strvcat(FilePath{NumFile}))
    %fprintf('\tOutput NODATA-Value:\t %f\n', nodata)
    %fprintf('\tSpatial Reference Mat:\t %s\n', ref)
    if ~isempty(type)
        %fprintf('\tGrid class type:\t %s\n', type)
    end
    
    start = 0;
    eval(['fid = fopen(''', strvcat(FilePath{NumFile}), ''',''r'');']);
    if fid==-1, error('Invalid file name!!'), end

    % ********************  Import HEADER  ****************************
    while start == 0;   % to be sure to read HEADER of file

        Read_header = fscanf(fid,'%s',1);
        % x- and y- llcenter are for instance used in ISATIS
        switch Read_header
            case {'ncols', 'NCOLS'}
                ncols = fscanf(fid,'%f',1);
            case {'nrows', 'NROWS'}
                nrows = fscanf(fid,'%f',1);
            case {'xllcorner', 'XLLCORNER', 'xllcenter'}
                xllcorner = fscanf(fid,'%f',1);
            case {'yllcorner', 'YLLCORNER', 'yllcenter'}
                yllcorner = fscanf(fid,'%f',1);
            case {'cellsize', 'CELLSIZE'}
                cellsize = fscanf(fid,'%f',1);
            case {'xdim'}
                xdim = fscanf(fid,'%f',1);
            case {'ydim'}
                ydim = fscanf(fid,'%f',1);                
            case {'NODATA_value', 'nodata_value', 'NODATA_VALUE'}
                NODATA_value = fscanf(fid,'%f',1);
                start = 1;
            otherwise
                error('Unrecognized %s in header',Read_header)
        end     %switch

    end         %while
    % *****************************************************************

    % ----------- pixel geometry -----------
    if exist('cellsize','var') % ==> the pixel is squared
        xdim = cellsize;
        ydim = cellsize;
    end
    % --------------------------------------  
    
    % -----------   Import RASTER   -----------
    import = fscanf(fid,'%f',[ncols nrows]);
    fclose(fid);
    % -----------------------------------------

    
    % +++++++++++++  Write No-Data Value  +++++++++++++++++++
    import(find(import == NODATA_value)) = nodata;
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++

    %mulsel = @(O,I,inv) [O '.' strrep(strvcat(FileName{NumFile}{:}(1:length(FileName{NumFile}{:})-4)),'-','_') ' = ' I inv ';'];
    mulsel = @(O,I,inv) [O '.' strrep(strvcat(FileName{NumFile}(1:length(FileName{NumFile})-4)),'-','_') ' = ' I inv ';'];
    % @@@@@@@@@@@@   Write RASTER   @@@@@@@@@@@@@@@
    if type == 's'                  % structure array of 2-D layers
        %eval(['Z.' strrep(strvcat(FileName{NumFile}{:}(1:length(FileName{NumFile}{:})-4)),'-','_') ' = import'';']);
        eval(mulsel('Z','import',''''))
    elseif type == 'd'              % double 3-D array stack
        Z(:,:,NumFile) = import';
    else                            % single 2-D layer
        Z = import';
    end
    % @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    
    % -----------------  Write HEADER  -------------------------
    if nargout == 2
        % Reference Matrix
        if ref == 'r'
            % ulc: Upper Left Centre point, that is the point in the centre of
            % the first upper-left pixel. This point is requested by
            % 'makerefmat.m' function.
            ulc_X = xllcorner + xdim /2;
            ulc_Y = yllcorner - ydim /2 + nrows*ydim;
            % This is the R spatial-referencing matrix used in Mapping Toolbox
            temp_R = makerefmat(ulc_X,ulc_Y,xdim,-ydim);
        % Header Matrix
        elseif ref == 'h'
            if ~exist('cellsize','var')
                warning('Unable to make canonical header')
                temp_R = [ncols; nrows; xllcorner; yllcorner; xdim; ydim; NODATA_value];
            else
                % This is the common header used in Arc/Info ASCII-Raster
                temp_R = [ncols; nrows; xllcorner; yllcorner; cellsize; NODATA_value];
            end
        end     %if
        % MultiSelection
        if type == 's'
            %eval(['R.' strrep(strvcat(FileName{NumFile}{:}(1:length(FileName{NumFile}{:})-4)),'-','_') ' = temp_R;']);
            eval(mulsel('R','temp_R',''))
        else
            R = temp_R;
        end
    end         %nargout
    % -----------------------------------------------------
    
    
end             %NumFile


%fprintf('...Completed!\n')

return

function [filename] = path2file(filepath)

fnd_b_slash = strfind(char(filepath), filesep);

if fnd_b_slash 
    STR = char(filepath);
    filename = cellstr(STR(fnd_b_slash(end)+1:end));
else
    filename = filepath;
end
return
