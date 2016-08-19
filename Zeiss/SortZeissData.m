function SortZeissData( filepath )
if ( nargin < 1 )
	path = input( 'Directory = ', 's' );
	fname = input( 'Filename = ', 's' );
else
	[ path, fname, ext ] = fileparts( filepath );
	fname = [ fname ext ];
end

%tzeiss = ReadZeissHeader( path, fname );
tzeiss = FindZeissHeader( path, fname );
if ( tzeiss.is_open )
	fprintf( 2, 'This file is sorted already: %s\\%s \n', path, fname );
	return;
end

fname = fname;
meta = ReadZeissMetaXML( fullfile( path, [ fname '_meta.xml' ] ) );
if ( ~meta.is_open )
	fprintf( 2, 'Opening META XML failed.\n' );
	return;
end
info = ReadZeissInfoXML( fullfile( path, [ fname '_info.xml' ] ) );
if ( ~meta.is_open )
	fprintf( 2, 'Opening INFO XML failed.\n' );
	return;
end

uscenes = [];
uscenes_StartS = [];

numZ = max( [ info.scenes.StartZ ] ) + 1;
numC = max( [ info.scenes.StartC ] ) + 1;
numT = max( [ info.scenes.StartT ] ) + 1;
numS = max( [ info.scenes.StartS ] ) + 1;
%maxB = max( [ info.scenes.StartB ] );
infomap = zeros( numZ, numC, numT, numS );
infomap_validated = false;
for i = 1:numel(info.scenes)
	scenei = info.scenes(i);
	%infomap( scenei.StartB+1, scenei.StartC+1, scenei.StartT+1, scenei.StartS+1 ) = i;
	infomap( scenei.StartZ+1, scenei.StartC+1, scenei.StartT+1, scenei.StartS+1 ) = i;

	if ~ismember( scenei.StartS, uscenes_StartS )
		uscenes = [ uscenes scenei ];
		uscenes_StartS(end+1) = scenei.StartS;
	end
end

czi = SortZeissCZI( fullfile( path, '..', [ fname '.czi' ] ) );
if czi.is_open
	meta.acqtime = czi.tdata_datenum(1) + czi.tdata_subsec(1)/24/3600 ;

	% XXX LKSCMT: Since it is unnecessarily difficult to extract the sorting order, it is hard-coded.
	if ( numel(czi.tdata_datenum) == numel(info.scenes) )
		infomap_validated = true;

		% Z - C - S - T
		i = 0;
		for tno = 1:numT
			for sno = 1:numS
				for cno = 1:numC
					for zno = 1:numZ
						i = i + 1;
						infomap( zno, cno, tno, sno ) = i;
					end
				end
			end
		end

	end
end


% Channels
nchannels = numel(meta.channels);
channels = repmat( struct, 1, nchannels );
ndigit = 1 + floor( log10( nchannels ) );
for i = 1:nchannels
	channels(i).name = meta.channels(i).name;
	channels(i).binX = meta.channels(i).binningX;
	channels(i).binY = meta.channels(i).binningY;
	channels(i).expT = meta.channels(i).exposureT;

	chstr = num2str(i);
	channels(i).name = [ 'CH' repmat( '0', 1, ndigit-length(chstr) ) chstr '.' channels(i).name ];
end

% Z-Stacks
nstacks = numel(meta.stacks);
stacks = repmat( struct, 1, nstacks );
ndigit = 1 + floor( log10( nstacks ) );
for i = 1:nstacks
	stacks(i).name = '';
	stacks(i).zpos = meta.stacks(i);

	if ( nstacks > 1 )
		stackstr = num2str(i);
		stacks(i).name = [ 'Z' repmat( '0', 1, ndigit-length(stackstr) ) stackstr ];
	end
end

% Scenes
[ nscenes, scenes, locations ] = ZeissScenes2Locations( meta.posArrays, uscenes );

destpath = cell( numel(locations), nchannels, nstacks );
for lcno = 1:numel(locations)
	lcpath = fullfile( path, locations(lcno).name );
	if ( exist( lcpath, 'file' ) <= 0 )
		mkdir( lcpath );
	end

	patches = locations(lcno).patches;
	% XXX LKSCMT: Assuming meta.px2umX == meta.px2umY. I don't like it...
	if ( WriteZeissHeader( lcpath, fname, meta.acqtime, meta.px2umX, meta.frame2ms, meta.nframes, channels, stacks, patches ) < 0 )
		fprintf( 2, 'Writing a header file failed.\n' );
		return;
	end
%	if czi.is_open
%		nchunk = numS * numZ * numC ;
%		lctI = ( numel(patches) * (lcno-1) ) + ( 1:nchunk:numel(czi.tdata_datenum) );
%		WriteZeissTData( lcpath, fname, czi.tdata_datenum(lctI) + czi.tdata_subsec(lctI)/24/3600 );		
%	end
	if infomap_validated
		lctI = reshape( infomap( 1, 1, :, 1 + numel(patches) * (lcno-1) ), 1, [] );
		WriteZeissTData( lcpath, fname, czi.tdata_datenum(lctI) + czi.tdata_subsec(lctI)/24/3600 );		
	end
	
	for chno = 1:nchannels
		chpath = fullfile( lcpath, channels(chno).name );
		if ( exist( chpath, 'file' ) <= 0 )
			mkdir( chpath );
		end
		for stno = 1:nstacks
			stpath = fullfile( chpath, stacks(stno).name );
			if ( exist( stpath, 'file' ) <= 0 )
				mkdir( stpath );
			end

			destpath{ lcno, chno, stno } = cell( 1, numel(patches) );
			for pcno = 1:numel(patches)
				destpath{ lcno, chno, stno }{pcno} = fullfile( stpath, patches(pcno).name );
				if ( exist( destpath{lcno, chno, stno }{pcno}, 'file' ) <= 0 )
					mkdir( destpath{lcno, chno, stno }{pcno} );
				end
			end
		end
	end
end

% Unix Timestamp
% ( Epoch time 1970/01/01 00:00:00 )
epoch = datenum( '1970-01-01 00:00:00' );
timestamp0 = 24 * 3600 * ( meta.acqtime - epoch );
% Java Timestamp
%timestamp0 = uint64(1000 * timestamp0);
timestamp0 = 1000 * timestamp0;

ndigit = 1 + floor( log10( meta.nframes ) );
files = dir( fullfile( path, [ fname '*' '_ORG.tif' ] ) );
for fno = 1:numel(files)
	tokens = regexp( files(fno).name, '(s\d+)*(t\d+)*(z\d+)*(c\d+)*_ORG.tif$', 'tokens' );
	lcno = 1;
	pcno = 1;
	chno = 1;
	stno = 1;

	if ( ~isempty(tokens) )
		stkn = [ 1 str2num(char(tokens{1}{1}(2:end))) ];
		lcno = scenes( stkn(end) ).lcno;
		pcno = scenes( stkn(end) ).pcno;
		ttkn = [ 1 str2num(char(tokens{1}{2}(2:end))) ];
		tstr = num2str( ttkn(end) );
		tstr = [ repmat( '0', 1, ndigit-length(tstr) ) tstr ];
		ztkn = [ 1 str2num(char(tokens{1}{3}(2:end))) ];
		stno = ztkn(end);
		ctkn = [ 1 str2num(char(tokens{1}{4}(2:end))) ];
		chno = ctkn(end);

		renameas = fullfile( destpath{ lcno, chno, stno }{pcno}, [ fname ' t' tstr '.tif' ] );
		java.io.File( fullfile( path, files(fno).name ) ).renameTo( java.io.File( renameas ) );

		timestamp = uint64( timestamp0 + meta.frame2ms * ( ttkn(end)-1 ) );
		if ( infomap_validated )
			cziIdx = infomap( ztkn(end), ctkn(end), ttkn(end), stkn(end) );
			timestamp = uint64( 1000 * ( 24 * 3600 * ( czi.tdata_datenum(cziIdx) - epoch ) + czi.tdata_subsec(cziIdx) ) );
		end
		java.io.File( renameas ).setLastModified( timestamp );
	end
end
return;

function tzeiss = FindZeissHeader( path, fname );
tzeiss = ReadZeissHeader(path, fname);
if ( tzeiss.is_open )
	fprintf( 2, 'This file is sorted already: %s\\%s \n', path, fname );
	return;
end

files = dir(path);
for i = 1:length(files)
	if ( ~strcmp(files(i).name, '.') && ~strcmp(files(i).name, '..') )
		tzeiss = FindZeissHeader( fullfile( path, files(i).name ), fname );
		if ( tzeiss.is_open )
			return;
		end
	end
end

return;
