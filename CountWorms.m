function CountWorms;
addpath( 'Zeiss' );

fclose('all');
path = input('Directory : ', 's');
fname = input('Filename : ', 's');

savpath = fullfile( path, [ fname '.KACnt.mat' ] );

global zeiss
if ( ~isempty(fname) )
	zeiss = OpenZeiss( fullfile( path, fname ), false );
else
	zeiss = OpenZeissDir( path, false );
end

w = zeiss.framew;
h = zeiss.frameh;
wh = max( w, h );

nframes = zeiss.nframes;
nchannels = zeiss.nchannels;
nstacks = zeiss.nstacks;

awbox = struct;
awbox.A = NaN;
awbox.CX = NaN;
awbox.CY = NaN;
awbox.P = NaN;
awbox.rect = zeros( 1, 4 );
awbox.worms = zeros( 0, 2 );
awbox.deads = zeros( 0, 2 );
awbox.notws = zeros( 0, 2 );
awbox.manual = false;

global cntdata
cntdata = struct;
cntdata.ROI = [ 1 1 ; w 1 ; w h ; 1 h ];
aframe = struct;
aframe.wboxes = repmat( awbox, 1, 0 );
cntdata.frames = repmat( aframe, 1, nframes );
cntdata.selected = 1:nframes;
if exist( savpath, 'file' ) > 0 
	% Load .KACnt.mat
	load( savpath, '-mat', 'cntdata' );
end


cmap = colormap(gray(256));

last_frameno = 0;
redraw_flag = true;
frameno = 1;
chno = 1;
stno = round(nstacks/2);
[ background, scaler ] = EstimateBS(chno, stno, frameno);
while 1
	frameno = min( frameno, nframes );
	frameno = max( 1, frameno );
	fprintf( 1, 'Frame %d of %d \n', frameno, nframes );
	fprintf( 1, '[ %d ] WormBoxes. \n', numel(cntdata.frames(frameno).wboxes) );

	if last_frameno ~= frameno
		frame = ReadZeiss( zeiss, chno, stno, frameno );
		[ temp, files ] = ReadZeiss( zeiss, chno, stno, frameno );
		[ temp, ffname ] = fileparts( files{1} );
		wpath = fullfile( path, [ ffname '.KACnt.png' ] );
		wframe = zeros( h, w );
		if exist( wpath, 'file' ) > 0 
			wframe = imread( wpath );
		end
		last_frameno = frameno;
		redraw_flag = true;
	end
	if redraw_flag
		uframe = max( 0.0, min( ( frame - background )/scaler, 1.0 ) );

		figure(1);
		clf(1);
		rgbframe = ones( h, w, 3 );
		xrgbframe = ones( wh, wh, 3 );
		rgbframe( :, :, 1 ) = max( wframe, uframe );
		rgbframe( :, :, 2 ) = ~wframe .* uframe;
		rgbframe( :, :, 3 ) = ~wframe .* uframe;
		xrgbframe( 1:h, 1:w, : ) = rgbframe;
		image( xrgbframe );
		axis off;
		axis image;
		zoom on;
		hold on;
		plot( cntdata.ROI([ 1:end 1 ], 1), cntdata.ROI([ 1:end 1 ], 2), ...
					'Color', [ 0 0 0 ], 'LineWidth', 3 );
		hold off;
		redraw_flag = false;
	end

	ans = input('[R]OI; Sa[V]e; [Q]uit ? ', 's');
	if strcmpi(ans, 'Q')
		disp('Quit');
		break;
	elseif strcmpi(ans, 'N')
		disp('Next');
		frameno = frameno + 1;
	elseif strcmpi(ans, 'P')
		disp('Previous');
		frameno = frameno - 1;
	elseif strcmpi(ans, 'G')
		disp('Go to');
		temp = round( input( sprintf( 'Frame No? [ 1 - %d ] ', nframes ) ) );
		if ( ~isempty(temp) )
			frameno = temp;
		end
	elseif strcmpi(ans, 'R')
		disp('Redefine ROI');
		ROI = zeros( 0, 2 );
		figure(1);
		hold on;
		hROI = plot( [ NaN ], [ NaN ], 'Color', [ 0 1 1 ], 'LineWidth', 3 );
		hold off;
		while 1
			figure(1);
			zoom on;
			ans = input( '[A]dd vertices; [U]se/[D]iscard this ROI ? ', 's' );
			if strcmpi( ans, 'U' ) 
				disp( 'The new ROI is defined.' );
				break;
			elseif strcmpi( ans, 'D' )
				disp( 'Discard this ROI' );
				break;
			elseif strcmpi( ans, 'A' )
				disp( 'Add more vertices' );
				figure(1);
				[ gx, gy ] = ginput(1);
				if ~isempty(gx)
					ROI = [ ROI ; gx gy ];
					set( hROI, 'XData', ROI([ 1:end 1 ], 1), ...
								'YData', ROI([ 1:end 1 ], 2) );
				end
			end
		end
		if strcmpi( ans, 'U' )
			fprintf( 1, 'Applying the new ROI...\n' );

			cntdata.ROI = ROI;
			for i = 1:nframes
				wboxes = cntdata.frames(i).wboxes;
				if ~isempty(wboxes)
					indi = find( ~inpolygon( [ wboxes.CX ], [ wboxes.CY ], ...
											ROI(:, 1), ROI(:, 2) ) );
					for k = reshape(indi, 1, [])
						wboxes(k).notws = [ wboxes(k).notws ; wboxes(k).worms ];
						wboxes(k).worms = zeros( 0, 2 );
					end
					cntdata.frames(i).wboxes = wboxes;
				end
			end
		end
		redraw_flag = true;
	elseif strcmpi(ans, 'V')
		disp('Save');
		asavpath = input( 'Save As : ', 's' );
		if ~isempty(asavpath)
			save( asavpath, 'cntdata' );
			fprintf( 1, 'CountData Saved As [%s].\n', asavpath );
		end
	elseif strcmpi(ans, 'BS')
		disp('Background/Scaler');
		fprintf( 1, 'Background: %.0f \n', background );
		fprintf( 1, 'Scaler: %.0f \n', scaler );
		background = [ background input( 'Background ? ' ) ];
		background = background(end);
		scaler = [ scaler input( 'Scaler ? ' ) ];
		scaler = scaler(end);
		redraw_flag = true;
	elseif strcmpi(ans, 'I')
		disp('Intensity Distribution');
		figure(2);
		clf(2);
		hold on;
		[ hN, hX ] = hist( reshape( frame, 1, [] ), 100 );
		plot( hX, hN, 'Color', [ 0 0.5 0 ], 'LineWidth', 2 );
		temp = axis;
		[ hN, hX ] = hist( reshape( dframep, 1, [] ), 100 );
		plot( hX, hN, 'Color', [ 0 0 1 ], 'LineWidth', 2 );
		[ hN, hX ] = hist( reshape( dframen, 1, [] ), 100 );
		plot( hX, hN, 'Color', [ 1 0 0 ], 'LineWidth', 2 );
		hold off;
		axis(temp);
		grid on;
	end
end
save( savpath, 'cntdata' );
fprintf( 1, 'CountData Saved.\n' );

CloseZeiss(zeiss);
return;

function [ background, scaler ] = EstimateBS( chno, stno, frameno );
global zeiss
framedata = sort( reshape( ReadZeiss( zeiss, chno, stno, frameno ), 1, [] ), 'ascend' );
background = framedata(ceil(0.01*end));
scaler = framedata(ceil(0.99*end)) - background;
%background = framedata(ceil(0.99*end));
%scaler = abs( framedata(ceil(0.01*end)) - background );
return;
