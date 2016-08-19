function CountWormsBatch( path, fname, redo );
addpath( 'Zeiss' );

fclose('all');
redo = true;
if nargin < 2
	path = input('Directory : ', 's');
	fname = input('Filename : ', 's');
elseif nargin < 3
	redo = false;
end

savdata = [];
savpath = fullfile( path, [ fname '.KACnt.mat' ] );
if exist( savpath, 'file' ) > 0 
	if ~redo
		fprintf( 2, 'Already Counted: %s \n', savpath );
		return;
	end
	% Load .KACnt.mat
	load( savpath, '-mat', 'cntdata' );
	savdata = cntdata;
end

global zeiss
if ( ~isempty(fname) )
	zeiss = OpenZeiss( fullfile( path, fname ), false );
else
	zeiss = OpenZeissDir( path, false );
end

w = zeiss.framew;
h = zeiss.frameh;
nframes = zeiss.nframes;
nchannels = zeiss.nchannels;
nstacks = zeiss.nstacks;


param = struct;
param.min_wormA = 300;	%Area
param.min_wormP = 100;	%Perimeter
param.wormR =  3;		%Radius
param.bgndR = 15;		%Background Radius
cntdata = struct;
cntdata.param = param;
cntdata.ROI = [ 1 1 ; w 1 ; w h ; 1 h ];
cntdata.frames = repmat( struct, 1, nframes );
cntdata.selected = 1:nframes;

if ~isempty(savdata)
	cntdata.ROI = savdata.ROI;
end

chno = 1;
stno = round(nstacks/2);
for frameno = 1:nframes
	frame = ReadZeiss( zeiss, chno, stno, frameno );
	%BW = MakeWormBWFrame( frame, param );
	[ cntdata.frames(frameno).wboxes, wBW ] = FindWormBoxes( frame, cntdata.param, cntdata.ROI );

	[ temp, files ] = ReadZeiss( zeiss, chno, stno, frameno );
	[ temp, ffname ] = fileparts( files{1} );
	filepath = fullfile( path, [ ffname '.KACnt.png' ] );
	imwrite( wBW, filepath, 'png', 'BitDepth', 1);
end

save( savpath, 'cntdata' );
CloseZeiss(zeiss);
return;
