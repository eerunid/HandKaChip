function zeiss = OpenZeissDir( path, readall );
zeiss = struct;
zeiss.is_open = false;

if ( nargin < 1 )
	path = input( 'Directory = ', 's' );
end

channel = struct;
channel.name = '';
%channel.binX = NaN;
%channel.binY = NaN;
%channel.expT = NaN;
channels = repmat( channel, 1, 0 );

subdirs = dir( path );
for i = 1:numel(subdirs)
	if ( subdirs(i).isdir && ~strcmp( subdirs(i).name, '.' ) && ~strcmp( subdirs(i).name, '..' ) )
		channel.name = subdirs(i).name;
		channels = [ channels channel ];
	end
end
if ( isempty(channels) )
	channel.name = '';
	channels = [ channels channel ];
end
nchannels = numel(channels);

files = dir( fullfile( path, channels(1).name, '*.tif' ) );
fnames = cell( 1, 0 );
fnamei = [];
lastfname = '';
lastfnamelen = 1;
for i = 1:numel(files)
	flag = true;
	for chno = 2:nchannels
		if ( exist( fullfile( path, channels(chno).name, files(i).name ), 'file' ) <= 0 )
			flag = false;
			break;
		end
	end

	if ( flag && ~strncmp( files(i).name, lastfname, lastfnamelen ) )
		tokens = regexp( files(i).name, '^(.*) t\d+.tif$', 'tokens' );
		%tokens = regexp( files(i).name, '^(.*)_ORG.tif$', 'tokens' );
		if ( ~isempty(tokens) )

			lastfname = char(tokens{1});
			lastfnamelen = length(lastfname);
			fnamei = [ fnamei i ];
			fnames = [ fnames lastfname ];
		end
	end
end
files = files(fnamei);
t = [ files(:).datenum ];
[ t, sortI ] = sort( t, 'ascend' );
files = files(sortI);
fnames = fnames(sortI);

if ( isempty(sortI) )
	return;
end


framew = NaN;
frameh = NaN;
delete = [];
for i = 1:numel(files)
	fullpathi = fullfile( path, channels(1).name, files(i).name );
	temp = dir( fullpathi );
	if ( temp.bytes > 0 )
		info = imfinfo( fullfile( path, channels(1).name, files(i).name ) );
		if ( ~isfinite(framew) || ~isfinite(frameh) )
			framew = info.Width;
			frameh = info.Height;
		end
		if ( framew == info.Width && frameh == info.Height )
			fprintf( 1, 'Reading [%s]...\n', files(i).name );
		else
			fprintf( 2, 'Discarding [%s] : ( %d, %d ) ~= ( %d, %d ) \n', ...
						files(i).name, info.Width, info.Height, framew, frameh );
			delete(end+1) = i;
		end
	end
end
fnames(delete) = [];

zeiss = struct;
zeiss.is_open = false;
zeiss.path = path;
zeiss.fname = '';
%zeiss.acqtime = NaN;
%zeiss.px2um = NaN;
%zeiss.frame2ms = NaN;

zeiss.nchannels = nchannels;
zeiss.channels = channels;

zeiss.nstacks = 1;
zeiss.stacks = repmat( struct, 1, 1 );
zeiss.stacks(1).name = '';
zeiss.stacks(1).zpos = NaN;

zeiss.npatches = 1;
zeiss.patches = repmat( struct, 1, 1 );
zeiss.patches(1).name = '';
zeiss.patches(1).l = 1;
zeiss.patches(1).r = framew;
zeiss.patches(1).b = 1;
zeiss.patches(1).t = frameh;

zeiss.framew = framew;
zeiss.frameh = frameh;

zeiss.nframes = 0;
zeiss.frames = [];
zeiss.tframes = [];
for i = 1:numel(fnames)
	ffiles = dir( fullfile( path, channels(1).name, [ fnames{i} '*.tif' ] ) );
	nffiles = numel(ffiles);

	fnomax = str2num( ffiles(end).name( (numel(fnames{i})+1):(end-4) ) );
	if ( fnomax ~= nffiles )
		fprintf( 2, '%s : Max. FrameNo[ %d ] ~= # of Files[ %d ] \n', fnames{i}, fnomax, nffiles );
	end

	zeiss.nframes = zeiss.nframes + nffiles ;
	zeiss.frames = [ zeiss.frames { ffiles.name } ];
	zeiss.tframes = [ zeiss.tframes [ ffiles.datenum ] ];
end
[ temp, sortI ] = sort( zeiss.tframes, 'ascend' );
zeiss.frames = zeiss.frames(sortI);
zeiss.tframes = ( zeiss.tframes(sortI) - zeiss.tframes(sortI(1)) ) * 24 * 60;

zeiss.is_open = true;
zeiss.bmask = uint8( zeros( 1, zeiss.nframes ) );
zeiss.bframes = cell( 1, zeiss.nframes );
if ( nargin >= 2 && readall )
	zeiss = LoadZeissBuffer( zeiss, 1:zeiss.nframes );
end
return;
