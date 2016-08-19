function zeiss = ReadZeissHeader( path, fname );
zeiss = struct;
zeiss.is_open = false;

%[ temp, fname, ext ] = fileparts( fname );
hdrpath = fullfile( path, [ fname '.dat' ] );
fid = fopen( hdrpath, 'r' );
if ( fid == -1 )
%	fprintf( 2, 'File open failed: %s \n', hdrpath );
	return;
end

zeiss.path = path;
%Filename
line = fgetl(fid);
line = fgetl(fid);	
if ( ~ischar(line) )
	fclose(fid);
	return;
end
zeiss.fname = line;

%Acquisition Time
line = fgetl(fid);
line = fgetl(fid);
if ( ~ischar(line) )
	fclose(fid);
	return;
end
zeiss.acqtime = str2num(line);
if ( numel(zeiss.acqtime) ~= 1 )
	fclose(fid);
	return;
end

% Length Per Pixel (um/px)
line = fgetl(fid);
line = fgetl(fid);
if ( ~ischar(line) )
	fclose(fid);
	return;
end
zeiss.px2um = str2num(line);
if ( numel(zeiss.px2um) ~= 1 )
	fclose(fid);
	return;
end

% Frame Rate (ms/frame)
line = fgetl(fid);
line = fgetl(fid);
if ( ~ischar(line) )
	fclose(fid);
	return;
end
zeiss.frame2ms = str2num(line);
if ( numel(zeiss.frame2ms) ~= 1 )
	fclose(fid);
	return;
end

% Number of Frames
line = fgetl(fid);
line = fgetl(fid);
if ( ~ischar(line) )
	fclose(fid);
	return;
end
zeiss.nframes = str2num(line);
zeiss.nframedigit = 1 + floor( log10( zeiss.nframes ) );
if ( numel(zeiss.nframes) ~= 1 || ~isfinite(zeiss.nframes) || zeiss.nframes < 0 )
	fclose(fid);
	return;
end

% Channels
%% Number of Channels
line = fgetl(fid);
line = fgetl(fid);
if ( ~ischar(line) )
	fclose(fid);
	return;
end
zeiss.nchannels = str2num(line);
if ( numel(zeiss.nchannels) ~= 1 || ~isfinite(zeiss.nchannels) || zeiss.nchannels < 0 )
	fclose(fid);
	return;
end
zeiss.channels = repmat( struct, 1, zeiss.nchannels );
for i = 1:zeiss.nchannels
	line = fgetl(fid);
	if ( ~ischar(line) )
		fclose(fid);
		return;
	end
	tokens = regexp( line, '^"(.*)" (\d+) (\d+) (\d+)$', 'tokens' );
	zeiss.channels(i).name = tokens{1}{1};
	zeiss.channels(i).binX = str2num( tokens{1}{2} );
	zeiss.channels(i).binY = str2num( tokens{1}{3} );
	zeiss.channels(i).expT = str2num( tokens{1}{4} );
end

% Z-Stacks
%% Number of Z-Stacks
line = fgetl(fid);
line = fgetl(fid);
if ( ~ischar(line) )
	fclose(fid);
	return;
end
zeiss.nstacks = str2num(line);
if ( numel(zeiss.nstacks) ~= 1 || ~isfinite(zeiss.nstacks) || zeiss.nstacks < 0 )
	fclose(fid);
	return;
end
zeiss.stacks = repmat( struct, 1, zeiss.nstacks );
ndigit = 1 + floor( log10( zeiss.nstacks ) );
for i = 1:zeiss.nstacks
	line = fgetl(fid);
	if ( ~ischar(line) )
		fclose(fid);
		return;
	end
	tokens = regexp( line, '^"(.*)" ([0-9.-+]+)$', 'tokens' );
	zeiss.stacks(i).name = tokens{1}{1};
	zeiss.stacks(i).zpos = str2num( tokens{1}{2} );
end

% Patches
%% Number of Patches
line = fgetl(fid);
line = fgetl(fid);
if ( ~ischar(line) )
	fclose(fid);
	return;
end
zeiss.npatches = str2num(line);
if ( numel(zeiss.npatches) ~= 1 || ~isfinite(zeiss.npatches) || zeiss.npatches < 0 )
	fclose(fid);
	return;
end
zeiss.patches = repmat( struct, 1, zeiss.npatches );
count = 0;
for i = 1:zeiss.npatches
	line = fgetl(fid);
	if ( ~ischar(line) )
		fclose(fid);
		return;
	end
	tokens = regexp( line, '^"(.*)" ([-0-9]+) ([-0-9]+) ([-0-9]+) ([-0-9]+)$', 'tokens' );
	zeiss.patches(i).name = tokens{1}{1};
	zeiss.patches(i).l = str2num( tokens{1}{2} );
	zeiss.patches(i).r = str2num( tokens{1}{3} );
	zeiss.patches(i).b = str2num( tokens{1}{4} );
	zeiss.patches(i).t = str2num( tokens{1}{5} );
	count = count + 1;
end
if ( count ~= zeiss.npatches )
	fprintf( 2, 'Invalid # Patches: [%d] ~= [%d] \n', zeiss.npatches, count );
	fclose(fid);
	return;
end
zeiss.framel = min( [ zeiss.patches(:).l ] ) ;
zeiss.frameb = min( [ zeiss.patches(:).b ] ) ;
for i = 1:zeiss.npatches
	zeiss.patches(i).l = zeiss.patches(i).l - zeiss.framel + 1;
	zeiss.patches(i).r = zeiss.patches(i).r - zeiss.framel + 1;
	zeiss.patches(i).b = zeiss.patches(i).b - zeiss.frameb + 1;
	zeiss.patches(i).t = zeiss.patches(i).t - zeiss.frameb + 1;
end

fclose(fid);

zeiss.framew = max( [ zeiss.patches(:).r ] );
zeiss.frameh = max( [ zeiss.patches(:).t ] );
zeiss.frames = cell( 1, 0 );
zeiss.is_open = true;

zeiss.bmask = uint8( zeros( 1, zeiss.nframes ) );
return;
