function [ tdata, flag ] = ReadZeissTData( path, fname, offset );
if ( nargin < 3 )
	offset = NaN;
end

tdata = [];
flag = false;

if exist( fullfile( path, [ fname '.tdata' ] ), 'file' ) > 0
	data = importdata( fullfile( path, [ fname '.tdata' ] ), '', 1 );
	acqtime = data.data;
	data = importdata( fullfile( path, [ fname '.tdata' ] ), '', 3 );
	tdata = data.data ;
	tdata = reshape( tdata / 60 / 24 + acqtime, 1, [] );	% LKSCMT: TData is in minutes as of 2014.12.23.
	flag = true;
else
	tzeiss = ReadZeissHeader( path, fname );
	if ( tzeiss.is_open )
		tdata = tzeiss.acqtime + ( 0:(tzeiss.nframes-1) ) * tzeiss.frame2ms / 24 / 3600 / 1000 ;

		ndigit = 1 + floor( log10( tzeiss.nframes ) );
		framestr = [ repmat( '0', 1, ndigit-1 ) '1' ];
		framefname = [ tzeiss.fname ' t' framestr '.tif' ];
		timestamp = FindTimestamp( path, framefname );
		if timestamp > 0
			if ~isfinite(offset)
				javaCal = java.util.Calendar.getInstance;
				offset = javaCal.get(javaCal.ZONE_OFFSET) + javaCal.get(javaCal.DST_OFFSET) ;
			end
			timestamp = timestamp - offset / 24 / 3600 / 1000 ;
			tdata = timestamp + ( 0:(tzeiss.nframes-1) ) * tzeiss.frame2ms / 24 / 3600 / 1000 ;
		end
	end
end
return;

function timestamp = FindTimestamp( path, fname );
fclose('all');
if isempty(path)
	return;
end
timestamp = NaN;
fullpath = fullfile( path, fname );
temp = dir( fullpath );
if ~isempty(temp)
	if ( temp.bytes > 0 )
		inf = imfinfo( fullpath );
		if ( isfield( inf, 'DigitalCamera' ) )
			dc = inf.DigitalCamera;
			if ( isfield( dc, 'DateTimeOriginal' ) && isfield( dc, 'SubsecTimeOriginal' ) )
				datetime = dc.DateTimeOriginal;
				subsec = dc.SubsecTimeOriginal;
				subsec = str2num(subsec) * 0.1^(length(subsec));
				timestamp = datenum(datetime, 'yyyy:mm:dd HH:MM:SS') + subsec / 24 / 3600 ;
				return;
			end
		end
	end	
end

files = dir(path);
for i = 1:length(files)
	if ( files(i).isdir )
		if ( ~strcmp(files(i).name, '.') && ~strcmp(files(i).name, '..') )
			timestamp = FindTimestamp( fullfile( path, files(i).name ), fname );
			if timestamp > 0
				return;
			end
		end
	end
end
fclose('all');
return;
