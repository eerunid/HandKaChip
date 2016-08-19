function zeiss = OpenZeiss( filepath, readall );
if ( nargin < 1 )
	path = input( 'Directory = ', 's' );
	fname = input( 'Filename = ', 's' );
else
	[ path, fname ] = fileparts( filepath );
end

zeiss = ReadZeissHeader( path, fname );
if ( ~zeiss.is_open )
	fprintf( 2, 'Opening the header file failed: %s \n', filepath );
	return;
end

zeiss.bgframe = zeros( zeiss.frameh, zeiss.framew, zeiss.nchannels, zeiss.nstacks );
for chno = 1:zeiss.nchannels
	chpath = fullfile( zeiss.path, zeiss.channels(chno).name );
	for stno = 1:zeiss.nstacks
		stpath = fullfile( chpath, zeiss.stacks(stno).name );
		%bgpath = fullfile( stpath, [ zeiss.fname '.max.tif' ] );
		bgpath = fullfile( stpath, [ zeiss.fname '.avg.tif' ] );
		if exist( bgpath, 'file' ) > 0
			frame = double( imread( bgpath ) );
			zeiss.bgframe( :, :, chno, stno ) = frame - max(frame(:));
		end
	end
end

tdata = ReadZeissTData( zeiss.path, zeiss.fname );
if isempty(tdata)
	zeiss.tframes = ( 0:(zeiss.nframes-1) ) * zeiss.frame2ms / 60000;
else
	zeiss.acqtime = tdata(1) ;
	zeiss.tframes = ( tdata - tdata(1) ) * 24 * 60 ;
end

zeiss.bmask = uint8( zeros( 1, zeiss.nframes ) );
zeiss.bframes = cell( 1, zeiss.nframes );
if ( nargin >= 2 && readall )
	zeiss = LoadZeissBuffer( zeiss, 1:zeiss.nframes );
end
return;
