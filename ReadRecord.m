%%-------------------------------------------------------
% [rec, info] = ReadRecord(recName, nOfFrames , startFrame)

% Input : 
%   recName - full path of folder with .tiff/.tif files or single .avi file, 
%                    or full path of .avi/.tif/.tiff file.
%                    Assuming gray scale image.
%   
%   nOfFrames  - [optional] read this number of frames. defualt = Inf
%   startFrame - [optional] first frame to read. default = 1
%
% Output : 
%   rec  - 3D matrix with the record, double
%   info - struct with three fields:
%          filename - parameters from filename
%          cam      - struct of the camera input parameters  
%          setup    - setup parameters as passed by the user to RecordFromCamera 
%
%%-------------------------------------------------------
function [rec, info] = ReadRecord( recName, nOfFrames , startFrame)
    
    %% Check input parameters
    if ~exist(recName,'file')
        error(['Record ''' recName ''' do not exist!'])
    end
    
    if ~exist('nOfFrames','var') || isempty(nOfFrames)
        nOfFrames = Inf ; % read all record
    end

    if ~exist('startFrame','var') || isempty(startFrame)
        startFrame = 1 ; % read all record
    end
    
    %% Read Record
    if exist(recName,'file') == 7 % it's a folder
        folderpath = recName;
        % find all .tiff or .tif files
        tiff_files = dir([folderpath, '\*.tiff']) ;
        avi_files  = dir([folderpath, '\*.avi']) ;
        if isempty(tiff_files) && isempty(avi_files)
            error(['There are no ''.tiff'' or ''.avi'' files in input folder ''' folderpath '''']);
        elseif ~isempty(tiff_files) && ~isempty(avi_files)
            error([ '''' folderpath ''' contains both ''.tiff'' or ''.avi'' files. It must contain only one the the file types ' ]);        
        elseif ~isempty(avi_files)
            if numel(avi_files) > 1
                error([ '"' folderpath '" must contain only one .avi file but contains ' num2str(numel(avi_files)) ' files.']);
            end
            [ rec , vH ]  = Avi2Matrix( fullfile(folderpath,avi_files.name) ,nOfFrames , startFrame);
            nBits = vH.BitsPerPixel;
            
        elseif ~isempty(tiff_files)
            [ rec , nBits ] = Tiff2Matrix(folderpath,nOfFrames,startFrame);              
        end
    else % it's a file    
        [~, ~ , ext] = fileparts(recName);
        if strcmp(ext,'.avi')
            [ rec , vH ] = Avi2Matrix( recName ,nOfFrames , startFrame);
            nBits = vH.BitsPerPixel;
        elseif strcmp(ext,'.tiff') || strcmp(ext,'.tif')
            t = Tiff(recName,'r');
            nBits = getTag(t,'BitsPerSample');
            rec = read(t);
            close(t);
            rec = rec(:,:,startFrame:startFrame+nOfFrames-1);                 
        else
            error(['Unsupported file type ' ext  ' . Supported types are .tif .tiff .avi '])
        end
    end

    rec = double(rec);
    if nBits == 16
        if all(mod(rec(1:400),16) == 0)
            rec = rec/16; % because Basler camera for some reason uses last 12 bits instead of first
        end
    end
    %% Read Info
    if nargout > 1
        info = GetRecordInfo(recName); 
        if exist('nBits','var')
            info.nBits = nBits;
            if info.nBits == 16
                info.nBits = 12; % basler camera has Mono12 or Mono8
            end
        end
    end
end


function  [rec , nBits] = Tiff2Matrix(filename, nOfFrames, startFrame)
    if ~exist('startFrame','var') || isempty(startFrame)
        startFrame = 1;
    end
        
    if ~exist('nOfFrames','var') || isempty(nOfFrames) 
        nOfFrames = Inf;
    end
    
    nBits = nan;
    
    if exist(filename,'file') == 7  % its a folder
        tiff_files = dir([filename filesep '\*.tiff']) ;
        if ~isinf(nOfFrames)
            if numel(tiff_files) - startFrame + 1 < nOfFrames
                error('File "%s" -> There is no enough frames (requested: %d , in file starting from frame %d: %d) ',filename,nOfFrames,startFrame,numel(tiff_files) - startFrame + 1);
            end
        else
            nOfFrames = numel(tiff_files) - startFrame + 1;
        end
                
        % get first image in order to find out the image size
        t = Tiff(fullfile(filename,tiff_files(1).name),'r');
        rec = nan(getTag(t,'ImageLength'),getTag(t,'ImageWidth'),nOfFrames);
        nBits = getTag(t,'BitsPerSample');
        close(t);

        % read all images
        for k = 1:nOfFrames
            t = Tiff(fullfile(filename,tiff_files(k+startFrame-1).name),'r');
            rec(:,:,k) = read(t);
            close(t);
        end
    elseif exist(filename,'file') == 2 && ( endsWith(filename,'.tiff') || endsWith(filename,'.tif') ) % its a file
        % TBD need to check, and implement more efficiently
%         rec = tiffreadVolume(filename,'PixelRegion',{[1 Inf],[1 Inf],[ startFrame (startFrame+nOfFrames-1)]); %available from MATLAB 2020b

        t = Tiff(filename,'r');
% %         offsets = getTag(t,'SubIFD')
% %         setDirectory(t,k);
% %         setSubDirectory(t,offsets(k));
        rec = read(t);
        nBits = getTag(t,'BitsPerSample');
        close(t);
  
        if ~isinf(nOfFrames)
            rec = rec(:,:,startFrame:(nOfFrames+startFrame-1));
        elseif startFrame > 1
            rec = rec(:,:,startFrame:end);
        end
    end
    
    
end


