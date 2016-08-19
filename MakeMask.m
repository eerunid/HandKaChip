function MakeMask;
addpath( 'Zeiss' );

fclose('all');
path = input('Directory : ', 's');
fname = input('Filename : ', 's');

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


global cntdata
cntpath = fullfile( path, [ fname '.KACnt.mat' ] );
if ~exist( cntpath, 'file' ) > 0 
	fprintf( 2, '[%s] Not Found.\n' );
	return;
end
% Load .KACnt.mat
load( cntpath, '-mat', 'cntdata' );

ROI = ones( h, w );
[ X, Y ] = meshgrid( 1:w, 1:h );
ROI(:) = 0;
in = inpolygon( X, Y, cntdata.ROI(:, 1), cntdata.ROI(:, 2) );
ROI(in) = 1;


% Hardcoded parameters
objsiz = 3;
chno = 1;	% Channel No
stno = 1;	% Stack No

[ X, Y ] = meshgrid( -objsiz:objsiz, -objsiz:objsiz );
kdisc = exp( -1.0 * ( X.^2 + Y.^2 )/(objsiz/3)^2 );
kroll = double( X.^2 + Y.^2 < 2.5^2 );

kdisc_th = 0.9 * sum(kdisc(:));


%LoadWormData();
tic;
for frameno = 1:nframes
	[ temp, files ] = ReadZeiss( zeiss, chno, stno, frameno );
	[ temp, ffname ] = fileparts( files{1} );
	wpath = fullfile( zeiss.path, [ ffname '.KACnt.png' ] );
	wframe = zeros( zeiss.frameh, zeiss.framew );
	if exist( wpath, 'file' ) > 0 
		wframe = imread( wpath );
	end

	CC = bwconncomp(wframe);
	zframe = zeros( size(wframe) );
	for i = 1:CC.NumObjects
		zframe( CC.PixelIdxList{i} ) = i;
	end

	wboxes = cntdata.frames(frameno).wboxes;
	for i = 1:numel(wboxes)
		if isempty(wboxes(i).worms) || ~isempty(wboxes(i).deads) || ~isempty(wboxes(i).notws)
%wboxes(i)
%ans = input( sprintf('FrameNo[%d/%d] WBox[%d/%d]', frameno, nframes, i, numel(wboxes) ) );
			indi = Rect2Indice( wboxes(i).rect );
			indx = setdiff( 1:numel(zframe), indi );
			zi = unique( zframe(indi) );
			for zz = reshape( zi, 1, [] )
				if isempty( find( zz == zframe(indx) ) )
					zframe(indi( find( zz == zframe(indi) ) )) = 0;
				end
			end
%image( zframe );
%axis image;
%ans = input('');
		end
	end

nw = min( 256, max(zframe(:)) );
figure(1);
clf(1);
image( zframe );
axis image;
colormap(hsv(nw));

	BW = zframe;
	BW1 = double( conv2( BW, kdisc, 'same' ) >= kdisc_th );
	BW1 = BW;
	BW2 = double( conv2( BW1, kroll, 'same' ) > 0 );

	mask = ROI;
	mask( find( ROI .* BW2 ) ) = 2;
	mask( find( ROI .* BW1 ) ) = 3;
	
	maskpath = fullfile( zeiss.path, [ ffname '.Mask.png' ] );
	imwrite( uint8(mask), maskpath, 'BitDepth', 8 );
	fprintf( 1, 'FrameNo[%d/%d] %.1f s...\n', frameno, nframes, toc );
end

CloseZeiss(zeiss);
return;

function ind = Rect2Indice(rect);
global zeiss
w = zeiss.framew;
h = zeiss.frameh;
ind = [];
for x = rect(1):rect(2)
	ind = [ ind h*(x-1) + ( rect(3):rect(4) ) ];
end
return;
