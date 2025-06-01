function info = GetRecordInfo(recordName)
    if ~exist(recordName,'file')
        error('Record "%s" does not exist',recordName);
    end
    
    if exist(recordName,'file') == 7 % it's a folder
        infoFilename = fullfile(recordName,'info.mat');
    else % its a file        
        infoFilename = [recordName(1:find(recordName=='.',1,'last')-1),'_info.mat'];
    end
    if exist(infoFilename,'file')
        info = load(infoFilename);
    end
    info.name = GetParamsFromFileName(recordName);
    
    %% Read records info
        %% Read Record
    if exist(recordName,'file') == 7 % it's a folder
        folderpath = recordName;
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
            info.fileType = '.avi';
            vH = VideoReader(fullfile(folderpath,avi_files(1).name) ); 
            info.nBits = vH.BitsPerPixel;           
        elseif ~isempty(tiff_files)
            tH = Tiff(fullfile(folderpath,tiff_files(1).name),'r');
            info.fileType = '.tiff';
            info.nBits = getTag(tH,'BitsPerSample'); 
            name_words = strsplit(tiff_files(1).name,'__');
            if numel(name_words) > 1 && ~isnan(str2double(name_words{2}))
                info.cameraSN = name_words{2};
            end    
        end
    else % it's a file    
        [~, ~ , ext] = fileparts(recordName);
        if strcmp(ext,'.avi')
            vH = VideoReader(recordName ); 
            info.fileType = '.avi';
            info.nBits = vH.BitsPerPixel;
        elseif strcmp(ext,'.tiff') || strcmp(ext,'.tif')
            tH = Tiff(recordName,'r');
            info.fileType = '.tiff';
            info.nBits = tH.BitsPerSample;            
            close(t);                   
        else
            error(['Unsupported file type ' ext  ' . Supported types are .tif .tiff .avi '])
        end
    end

    if info.nBits == 16
        info.nBits = 12;
    end
    
    
    
