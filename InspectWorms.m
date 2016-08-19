function InspectWorms;
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
nframes = zeiss.nframes;
nchannels = zeiss.nchannels;
nstacks = zeiss.nstacks;

awbox = struct;
awbox.A = NaN;
awbox.CX = NaN;
awbox.CY = NaN;
awbox.P = NaN;
awbox.D = NaN;
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

global frame
frame = struct;
frame.curr_frameno = 0;
frame.curr_frame = zeros( w, h );
frame.curr_wframe = zeros( w, h );
frame.chno = 1;
frame.stno = round(nstacks/2);

LoadWormData(true);
redraw_flag = true;
frameno = 1;
[ frame.background, frame.scaler ] = EstimateBS( frame.chno, frame.stno, frameno );
while 1
	frameno = min( frameno, nframes );
	frameno = max( 1, frameno );

	redraw_flag = LoadFrame(frameno) || redraw_flag;
	if ( redraw_flag )
		DrawFrame();
		[ hCurrA, hCurrD ] = PlotWormData();
		redraw_flag = false;
	end
	wdata = [];

	ans = input('Filter By [A]rea/[D]arkness; Sa[V]e; [Q]uit ? ', 's');
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
	elseif strcmpi(ans, 'S')
		disp('Select Frames');
		fprintf( 1, 'Selected : %s \n', num2str(cntdata.selected) );
		fprintf( 2, 'Unselected : %s \n', num2str(setdiff( 1:nframes, cntdata.selected)) );
		ans = input( 'Select : ', 's' );
		flist = intersect( 1:nframes, eval( [ '[ ' ans ' ]' ] ) );
		cntdata.selected = reshape( union( cntdata.selected, flist ), 1, [] );
		LoadWormData(true);
		redraw_flag = true;
	elseif strcmpi(ans, 'U')
		disp('Unselect Frames');
		fprintf( 1, 'Selected : %s \n', num2str(cntdata.selected) );
		fprintf( 2, 'Unselected : %s \n', num2str(setdiff( 1:nframes, cntdata.selected)) );
		ans = input( 'Unselect : ', 's' );
		flist = intersect( 1:nframes, eval( [ '[ ' ans ' ]' ] ) );
		cntdata.selected = setdiff( cntdata.selected, flist );		
		LoadWormData(true);
		redraw_flag = true;
	elseif strcmpi(ans, 'W')
		disp('Add A WormBox');
		figure(1);
		k = waitforbuttonpress;
		point1 = get(gca, 'CurrentPoint');    % button down detected
		finalRect = rbbox;                   % return figure units
		point2 = get(gca, 'CurrentPoint');    % button up detected
		gx = sort( [ point1( 1, 1 ) point2( 1, 1 ) ], 'ascend' );
		gy = sort( [ point1( 1, 2 ) point2( 1, 2 ) ], 'ascend' );
		if ( gx(1) < gx(2) || gy(1) < gy(2) )
			rect = [ gx(1) gx(2) gy(1) gy(2) ];
			fprintf( 1, 'Add a new WormBox here ? \n' );
			fprintf( 1, 'Click the Worm. \n' );
			[ gx, gy ] = ginput(1);
			if ( numel(gx) > 0 ...
				&& rect(1) <= gx(1) && gx(1) <= rect(2) ...
				&& rect(3) <= gy(1) && gy(1) <= rect(4) )
				newwbox = awbox;
				newwbox.CX = gx(1);
				newwbox.CY = gy(1);
				newwbox.rect = rect;
				newwbox.worms = [ gx(1) gy(1) ];
				newwbox.manual = true;
				cntdata.frames(frameno).wboxes(end+1) = newwbox;
				redraw_flag = true;
			end
		end
	elseif strcmpi(ans, 'A')
		disp('Filter By Area');
		disp('Select ROI Polygon');
		figure(3);
		[ xv, yv ] = SelectPoly();
		if ~isempty(xv)
			wdata = LoadWormData(false);
			wdata = wdata( find( inpolygon( [ wdata.A ], [ wdata.frameno ], xv, yv ) ) );
			wdata( find( [ wdata.manual ] ) ) = [];
			fprintf( 1, '%d Worms Selected \n', numel(wdata) );
			ans = input( '[A]uto Fix/Manual [C]heck ? ', 's' );
			if strcmpi( ans, 'A' )
				disp( 'Auto Fix' );
				disp( 'Click for the Standard Size for a Single Worm' );
				figure(3);
				[ gx, gy ] = ginput(1);
				if numel(gx) > 0 && gx(1) > 0
					stdA = gx(1);
					wdata = wdata( find( [ wdata.ndeads ] + [ wdata.nnotws ] <= 0 ) );
					for wno = 1:numel(wdata)
						wfno = wdata(wno).frameno;
						wwno = wdata(wno).wboxno;
						wbox = cntdata.frames(wfno).wboxes(wwno);
						if ( wdata(wno).A < stdA )
							wbox.notws = [ wbox.notws ; wbox.worms ];
							wbox.worms = zeros( 0, 2 );
						else
							wn = round( wbox.A/stdA );
							wx = linspace( wbox.rect(1), wbox.rect(2), wn+2 );
							wy = linspace( wbox.rect(3), wbox.rect(4), wn+2 );
							wbox.worms = [ wx(2:end-1)' wy(2:end-1)' ];
						end
						cntdata.frames(wfno).wboxes(wwno) = wbox;
					end
					LoadWormData(true);
					redraw_flag = true;
				end
				wdata = [];
			elseif strcmpi( ans, 'C' )
				disp( 'Manual Check' );
			else
				wdata = [];
			end
		end
	elseif strcmpi(ans, 'D')
		disp('Filter By Darkness');
		disp('Select ROI Polygon');
		figure(4);
		[ xv, yv ] = SelectPoly();
		if ~isempty(xv)
			wdata = LoadWormData(false);
			wdata = wdata( find( inpolygon( [ wdata.D ], [ wdata.frameno ], xv, yv ) ) );
			wdata( find( [ wdata.manual ] ) ) = [];
			fprintf( 1, '%d Worms Selected \n', numel(wdata) );
			ans = input( '[A]uto Fix/Manual [C]heck ? ', 's' );
			if strcmpi( ans, 'A' )
				disp( 'Auto Fix' );
				disp( 'Selected Will be Discarded.' );
				wdata = wdata( find( [ wdata.ndeads ] + [ wdata.nnotws ] <= 0 ) );
				for wno = 1:numel(wdata)
					wfno = wdata(wno).frameno;
					wwno = wdata(wno).wboxno;
					wbox = cntdata.frames(wfno).wboxes(wwno);
					wbox.notws = [ wbox.notws ; wbox.worms ];
					wbox.worms = zeros( 0, 2 );
					cntdata.frames(wfno).wboxes(wwno) = wbox;
				end
				wdata = [];
				LoadWormData(true);
				redraw_flag = true;
			elseif strcmpi( ans, 'C' )
				disp( 'Manual Check' );
			else
				wdata = [];
			end
		end
	elseif strcmpi(ans, 'R')
		disp('Filter By ROI');
		disp(' - Select Frames');
		ans = input( 'Select : ', 's' );
		flist = intersect( 1:nframes, eval( [ '[ ' ans ' ]' ] ) );
		disp(' - Select ROI');

		figure(1);
		k = waitforbuttonpress;
		point1 = get(gca, 'CurrentPoint');    % button down detected
		finalRect = rbbox;                   % return figure units
		point2 = get(gca, 'CurrentPoint');    % button up detected
		gx = sort( [ point1( 1, 1 ) point2( 1, 1 ) ], 'ascend' );
		gy = sort( [ point1( 1, 2 ) point2( 1, 2 ) ], 'ascend' );
		if ( gx(1) < gx(2) && gy(1) < gy(2) )
			xv = gx( [ 1 2 2 1 ] );
			yv = gy( [ 1 1 2 2 ] );
			wdata = LoadWormData(false);
			wdata = wdata( find( inpolygon( [ wdata.CX ], [ wdata.CY ], xv, yv ) ) );
			if ~isempty(flist)
				wflag = ismember( [ wdata.frameno ], flist );
				wdata = wdata( find(wflag) );
			end
			fprintf( 1, '%d Worms Selected \n', numel(wdata) );
			ans = input( '[F]ilter them out first ? ', 's' );
			if strcmpi( ans, 'F' )
				ind = find( ~[ wdata.manual ] );
				for wno = reshape( ind, 1, [] )
					wfno = wdata(wno).frameno;
					wwno = wdata(wno).wboxno;
					wbox = cntdata.frames(wfno).wboxes(wwno);
					wbox.notws = [ wbox.notws ; wbox.worms ];
					wbox.worms = zeros( 0, 2 );
					cntdata.frames(wfno).wboxes(wwno) = wbox;
				end
			end
		end
	elseif strcmpi(ans, 'K')
		disp('Manual Pick');
		figure(1);
		[ gx, gy ] = ginput(1);
		if ~isempty(gx)
			wdata = LoadWormData(false);
			wrect = cat( 1, wdata.rect );
			wtest = ( [ wdata.frameno ]' == frameno ) ...
					.* ( (wrect(:, 1)-gx(1)).*(wrect(:, 2)-gx(1)) <= 0 ) ...
					.* ( (wrect(:, 3)-gy(1)).*(wrect(:, 4)-gy(1)) <= 0 );
			wdata = wdata( find( wtest ) );
			fprintf( 1, '%d Worms Selected \n', numel(wdata) );
		end
	elseif strcmpi(ans, 'V')
		disp('Save');
		save( savpath, 'cntdata' );
		fprintf( 1, 'CountData Saved.\n' );
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
		figure(5);
		clf(5);
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
	
	if ~isempty(wdata)
		CheckWorms( wdata, [ hCurrA hCurrD ] );
		LoadWormData(true);
		redraw_flag = true;
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

function h = DeleteHandle(h);
for h = reshape( h, 1, [] )
	if ishandle(h)
		delete(h);
	end
end
h = [];
return;

function MarkWorms( wboxes );
if isempty(wboxes)
	return;
end
global frame
frame.hWorms = DeleteHandle(frame.hWorms);

currFig = gcf;

worms = cat( 1, wboxes.worms );
deads = cat( 1, wboxes.deads );
notws = cat( 1, wboxes.notws );

figure(1);
hold on;
h1 = plot( worms(:, 1), worms(:, 2), 'o', 'Color', [ 0 1 0 ], 'LineWidth', 2, 'MarkerSize', 5 );
h2 = plot( deads(:, 1), deads(:, 2), 'o', 'Color', [ 0 1 1 ], 'LineWidth', 2, 'MarkerSize', 10 );
h3 = plot( notws(:, 1), notws(:, 2), 'x', 'Color', [ 0 1 1 ], 'LineWidth', 2, 'MarkerSize', 10 );
hold off;
frame.hWorms = [ h1' h2' h3' ];

figure( currFig );
return;

function MarkBox( wbox, rects );
global frame
frame.hBox = DeleteHandle(frame.hBox);

currFig = gcf;

figure(1);
rectX = [ wbox.rect(1) wbox.rect(2) wbox.rect(2) wbox.rect(1) wbox.rect(1) ];
rectY = [ wbox.rect(3) wbox.rect(3) wbox.rect(4) wbox.rect(4) wbox.rect(3) ];
hold on;
h = plot( rectX, rectY, 'Color', [ 1 1 0 ], 'LineWidth', 2 );
hold off;
frame.hBox(end+1) = h;

figure(2);
x0 = wbox.rect(1)-1;
y0 = wbox.rect(3)-1;
hold on;
% Other Rects
rects( find( wbox.rect(2) < rects(:, 1) ), : ) = [];
rects( find( rects(:, 2) < wbox.rect(1) ), : ) = [];
rects( find( wbox.rect(4) < rects(:, 3) ), : ) = [];
rects( find( rects(:, 4) < wbox.rect(3) ), : ) = [];
for i = 1:size(rects, 1)
	rectXi = [ rects(i, 1) rects(i, 2) rects(i, 2) rects(i, 1) rects(i, 1) ]-x0;
	rectYi = [ rects(i, 3) rects(i, 3) rects(i, 4) rects(i, 4) rects(i, 3) ]-y0;
	h = plot( rectXi, rectYi, 'Color', [ 0 0 0 ], 'LineWidth', 2 );
	frame.hBox(end+1) = h;
end

h1 = plot( wbox.worms(:, 1)-x0, wbox.worms(:, 2)-y0, 'o', ...
			'Color', [ 0 1 0 ], 'LineWidth', 2, 'MarkerSize', 10 );
h2 = plot( wbox.deads(:, 1)-x0, wbox.deads(:, 2)-y0, 'o', ...
			'Color', [ 0 1 1 ], 'LineWidth', 2, 'MarkerSize', 10 );
h3 = plot( wbox.notws(:, 1)-x0, wbox.notws(:, 2)-y0, 'x', ...
			'Color', [ 0 1 1 ], 'LineWidth', 2, 'MarkerSize', 10 );
hold off;
frame.hBox = [ frame.hBox h1 h2 h3 ];

figure( currFig );
return;


function reload_flag = LoadFrame( frameno )
global zeiss
global frame
reload_flag = ( frame.curr_frameno ~= frameno );
if reload_flag
	zframe = ReadZeiss( zeiss, frame.chno, frame.stno, frameno );
	[ temp, files ] = ReadZeiss( zeiss, frame.chno, frame.stno, frameno );
	[ temp, ffname ] = fileparts( files{1} );
	wpath = fullfile( zeiss.path, [ ffname '.KACnt.png' ] );
	wframe = zeros( size(zframe) );
	if exist( wpath, 'file' ) > 0 
		wframe = imread( wpath );
	end

	frame.w = zeiss.framew;
	frame.h = zeiss.frameh;
	frame.wh = max( frame.w, frame.h );
	frame.curr_frameno = frameno;
	frame.curr_zframe = zframe;
	frame.curr_wframe = wframe;
	frame.curr_uframe = zeros( size(zframe) );
end
return;

function DrawFrame( )
global cntdata
global frame

wframe = frame.curr_wframe;
uframe = max( 0.0, min( ( frame.curr_zframe - frame.background )/frame.scaler, 1.0 ) );
frame.curr_uframe = uframe;
frame.hWorms = [];
frame.hBox = [];

frameno = frame.curr_frameno;
is_selected = ismember(frameno, cntdata.selected);
fprintf( 1+~is_selected, 'Frame %d of %d \n', frameno, numel(cntdata.frames) );
wboxes = cntdata.frames(frameno).wboxes;
fprintf( 1+~is_selected, '[ %d ] WormBoxes. \n', numel(wboxes) );
fprintf( 1+~is_selected, '[ %d ] Worms. \n', size( cat(1, wboxes.worms), 1) );

figure(1);
clf(1);
rgbframe = ones( frame.h, frame.w, 3 );
xrgbframe = ones( frame.wh, frame.wh, 3 );
rgbframe( :, :, 1 ) = max( wframe, uframe );
rgbframe( :, :, 2 ) = ~wframe .* uframe;
rgbframe( :, :, 3 ) = ~wframe .* uframe;
xrgbframe( 1:frame.h, 1:frame.w, : ) = rgbframe;
image( xrgbframe );
axis off;
axis image;
zoom on;
hold on;
plot( cntdata.ROI([ 1:end 1 ], 1), cntdata.ROI([ 1:end 1 ], 2), ...
					'Color', [ 0 0 0 ], 'LineWidth', 3 );
hold off;
MarkWorms( cntdata.frames(frame.curr_frameno).wboxes );
return;

function wdata = LoadWormData(reload_flag);
global curr_wdata
if ~reload_flag
	wdata = curr_wdata;
	return;
end

global cntdata

aworm = struct;
aworm.frameno = 0;
aworm.wboxno = 0;
aworm.A = NaN;
aworm.D = NaN;
aworm.CX = NaN;
aworm.CY = NaN;
aworm.rect = zeros( 1, 4 );
aworm.nworms = 0;
aworm.ndeads = 0;
aworm.nnotws = 0;
aworm.manual = false;

wdata = repmat( aworm, 1, 0 );
for i = 1:numel(cntdata.frames)
	wboxesi = cntdata.frames(i).wboxes;
	for k = 1:numel(wboxesi)
		wboxk = wboxesi(k);
		wdata(end+1) = aworm;
		wdata(end).frameno = i;
		wdata(end).wboxno = k;
		wdata(end).A = wboxk.A;
		wdata(end).D = wboxk.D;
		wdata(end).CX = wboxk.CX;
		wdata(end).CY = wboxk.CY;
		wdata(end).rect = wboxk.rect;
		wdata(end).nworms = size( wboxk.worms, 1 );
		wdata(end).ndeads = size( wboxk.deads, 1 );
		wdata(end).nnotws = size( wboxk.notws, 1 );
		wdata(end).manual = wboxk.manual;
	end
end

% Only care for selected frames
wdata = wdata( find( ismember( [ wdata.frameno ], cntdata.selected ) ) );

curr_wdata = wdata;
return;


function [ hArea, hDarkness ] = PlotWormData();
global frame
global cntdata
wdata = LoadWormData(false);

%wdata = wdata( find( ismember( [ wdata.frameno ], cntdata.selected ) ) );

% Dead Worms and Not Worms
indN = intersect( find( [ wdata.ndeads ] + [ wdata.nnotws ] > 0 ), ...
					find( [ wdata.nworms ] <= 0 ) );
indN = union( find( [ wdata.ndeads ] + [ wdata.nnotws ] > 0 ), ...
				find( [ wdata.nworms ] <= 0 ) );
% Live Worms
indW = intersect( find( [ wdata.ndeads ] + [ wdata.nnotws ] <= 0 ), ...
					find( [ wdata.nworms ] > 0 ) );
indC = indW( find( [ wdata(indW).frameno ] == frame.curr_frameno ) );
indM = indW( find( [ wdata(indW).manual ] ) );

fdata = [ wdata.frameno ];
%adata = [ wdata.A ]./[ wdata.nworms ];
adata = [ wdata.A ]./max( 1, [ wdata.nworms ] );
ddata = [ wdata.D ];

%adata( find( [ wdata.nworms ] <= 0 ) ) = NaN;
adata( intersect( find( [ wdata.nworms ] > 1 ), ...
					find( ~[ wdata.manual ] ) ) ) = NaN;

	
figure(3);
clf(3);
hAx31 = subplot( 2, 1, 1 );
hold on;
[ hY, hX ] = hist( adata(indW), 100 );
bar( hX, hY, 'FaceColor', 0.5 * [ 1 1 1 ], 'EdgeColor', [ 1 1 1 ] );
hY = hist( adata(indM), hX );
bar( hX, hY, 'FaceColor', [ 0 0 0 ], 'EdgeColor', [ 1 1 1 ] );
hY = hist( adata(indC), hX );
bar( hX, hY, 'FaceColor', [ 1 0 0 ], 'EdgeColor', [ 1 1 1 ] );
hold off;
xlabel( 'Area (px)' );
hAx32 = subplot( 2, 1, 2 );
hold on;
plot( adata(indW), fdata(indW), 'x', 'Color', 0.5 * [ 1 1 1 ], 'MarkerSize', 10 );
plot( adata(indN), fdata(indN), 'x', 'Color', [ 0 1 1 ], 'LineWidth', 2, 'MarkerSize', 10 );
plot( adata(indM), fdata(indM), 'x', 'Color', [ 0 0 0 ], 'LineWidth', 2, 'MarkerSize', 10 );
plot( adata(indC), fdata(indC), 'x', 'Color', [ 1 0 0 ], 'LineWidth', 2, 'MarkerSize', 10 );
hArea = plot( [ NaN ], [ NaN ], 'o', 'Color', [ 0 1 0 ], 'LineWidth', 2, 'MarkerSize', 10 );
hold off;
xlabel( 'Area (px)' );
ylabel( 'Frame' );
linkaxes( [ hAx31 hAx32 ], 'x' );
figure(4);
clf(4);
hAx41 = subplot( 2, 1, 1 );
hold on;
[ hY, hX ] = hist( ddata(indW), 100 );
bar( hX, hY, 'FaceColor', 0.5 * [ 1 1 1 ], 'EdgeColor', [ 1 1 1 ] );
hY = hist( ddata(indM), hX );
bar( hX, hY, 'FaceColor', [ 0 0 0 ], 'EdgeColor', [ 1 1 1 ] );
hY = hist( ddata(indC), hX );
bar( hX, hY, 'FaceColor', [ 1 0 0 ], 'EdgeColor', [ 1 1 1 ] );
hold off;
xlabel( 'Worm-Brightness' );
hAx42 = subplot( 2, 1, 2 );
hold on;
plot( ddata(indW), fdata(indW), 'x', 'Color', 0.5 * [ 1 1 1 ], 'MarkerSize', 10 );
plot( ddata(indN), fdata(indN), 'x', 'Color', [ 0 1 1 ], 'LineWidth', 2, 'MarkerSize', 10 );
plot( ddata(indM), fdata(indM), 'x', 'Color', [ 0 0 0 ], 'LineWidth', 2, 'MarkerSize', 10 );
plot( ddata(indC), fdata(indC), 'x', 'Color', [ 1 0 0 ], 'LineWidth', 2, 'MarkerSize', 10 );
hDarkness = plot( [ NaN ], [ NaN ], 'o', 'Color', [ 0 1 0 ], 'LineWidth', 2, 'MarkerSize', 10 );
hold off;
xlabel( 'Worm-Brightness' );
ylabel( 'Frame' );
linkaxes( [ hAx41 hAx42 ], 'x' );
return;



function [ wormI, deadI, notwI ] = SelectWorms( wbox );
wormI = [];
deadI = [];
notwI = [];

figure(2);
k = waitforbuttonpress;
point1 = get(gca, 'CurrentPoint');    % button down detected
finalRect = rbbox;                   % return figure units
point2 = get(gca, 'CurrentPoint');    % button up detected
gx = sort( [ point1( 1, 1 ) point2( 1, 1 ) ], 'ascend' );
gy = sort( [ point1( 1, 2 ) point2( 1, 2 ) ], 'ascend' );
if ( gx(1) >= gx(2) || gy(1) >= gy(2) )
	return;
end
gx = gx + wbox.rect(1)-1;
gy = gy + wbox.rect(3)-1;
xv = gx( [ 1 2 2 1 ] );
yv = gy( [ 1 1 2 2 ] );
wormI = find( inpolygon( wbox.worms(:, 1), wbox.worms(:, 2), xv, yv ) );
deadI = find( inpolygon( wbox.deads(:, 1), wbox.deads(:, 2), xv, yv ) );
notwI = find( inpolygon( wbox.notws(:, 1), wbox.notws(:, 2), xv, yv ) );
return;

function [ nclicks, wormI, deadI, notwI, newpos ] = ClickAWorm( wbox );
nclicks = 0;
wormI = [];
deadI = [];
notwI = [];
newpos = zeros( 0, 2 );

figure(2);
[ gx, gy ] = ginput(1);
if ( numel(gx) <= 0 )
	return;
end
nclicks = 1;
gx = gx(1) + wbox.rect(1)-1;
gy = gy(1) + wbox.rect(3)-1;

objs = [ wbox.worms ; wbox.deads ; wbox.notws ];
[ dist2, minI ] = min( ( objs(:, 1) - gx ).^2 + ( objs(:, 2) - gy ).^2 );

global cntdata
if ( dist2 <= cntdata.param.wormR^2 )
	nworms = size(wbox.worms, 1);
	ndeads = size(wbox.deads, 1);
	if ( minI <= nworms )
		wormI = minI;
	elseif ( minI <= nworms + ndeads )
		deadI = minI - nworms;
	else
		notwI = minI - (nworms+ndeads);
	end
else
	newpos = [ gx gy ];
end

return;

function CheckWorms( wdata, hCurrent );
global cntdata
global frame

hA = hCurrent(1);
hD = hCurrent(2);
nwdata = numel(wdata);
i = 1;
while 1
	i = min( i, nwdata );
	i = max( 1, i );

	wfno = wdata(i).frameno;
	wwno = wdata(i).wboxno;
	if LoadFrame( wfno )
		DrawFrame();
	end
	wboxi = cntdata.frames( wfno ).wboxes(wwno);
	li = floor( wboxi.rect(1) );
	ri =  ceil( wboxi.rect(2) );
	bi = floor( wboxi.rect(3) );
	ti =  ceil( wboxi.rect(4) );
	nwormsi = size(wboxi.worms, 1);
	ndeadsi = size(wboxi.deads, 1);
	nnotwsi = size(wboxi.notws, 1);
	fprintf( 1, 'WormBox[ %d / %d ] ( %.0f-%.0f, %.0f-%.0f ) at Frame %d\n', ...
						i, nwdata, wboxi.rect, wfno );
	fprintf( 1, '	A[%.1f] P[%.1f] D[%.1f] \n', wboxi.A, wboxi.P, wboxi.D );
	for k = 1:nwormsi
		fprintf( 1, 'Worm[%d] ( %.0f, %.0f ) \n', k, wboxi.worms(k, :) );
	end
	for k = 1:ndeadsi
		fprintf( 2, 'Dead[%d] ( %.0f, %.0f ) \n', k, wboxi.deads(k, :) );
	end
	for k = 1:nnotwsi
		fprintf( 2, '"Not"Worms[%d] ( %.0f, %.0f ) \n', k, wboxi.notws(k, :) );
	end
	
	figure(2);
	clf(2);
	image( uint8( 255 * frame.curr_uframe( bi:ti, li:ri ) ) );
	colormap(gray(256));
	axis off;
	axis image;
	zoom on;
	brdata = cat( 1, cntdata.frames(wfno).wboxes.rect );
	MarkBox( wboxi, brdata );

	set( hA, 'XData', wboxi.A, 'YData', wfno );
	set( hD, 'XData', wboxi.D, 'YData', wfno );
	
	ans = input('[A]dd/[D]elete Worms/NotWorms; Pick [W]orm/D[E]ad/[X]Not Worms; [Q]uit ? ', 's');
	if strcmpi(ans, 'Q')
		disp('Quit');
		break;
	elseif strcmpi(ans, 'N')
		disp('Next');
		i = i + 1;
	elseif strcmpi(ans, 'P')
		disp('Previous');
		i = i - 1;
	elseif strcmpi(ans, 'G')
		disp('Go to');
		temp = round( input( sprintf( 'WormBox No ? [ 1 - %d ] ', nwdata ) ) );
		if ( ~isempty(temp) )
			i = temp;
		end
	elseif strcmpi(ans, 'M')
		disp('Next');
		fprintf( 2, 'This worm is manually approved.\n' );
		wboxi.manual = true;
		i = i + 1;
	elseif strcmpi(ans, 'A')
		disp('Add Worms/NotWorms');
		while 1
			figure(2);
			[ gx, gy ] = ginput(1);
			if ( numel(gx) <= 0 )
				break;
			end
			gx = gx(1) + wboxi.rect(1)-1;
			gy = gy(1) + wboxi.rect(3)-1;
			wboxi.worms( end+1, : ) = [ gx gy ];
			wboxi.manual = true;
			MarkBox( wboxi, brdata );
		end
	elseif strcmpi(ans, 'D')
		disp('Delete Worms/NotWorms');
		[ wormI, deadI, notwI ] = SelectWorms( wboxi );
		wboxi.worms(wormI, :) = [];
		wboxi.deads(deadI, :) = [];
		wboxi.notws(notwI, :) = [];
		wboxi.manual = true;
		MarkBox( wboxi, brdata );
	elseif strcmpi(ans, 'W')
		disp('Pick Worms');
		while 1
			[ nclicks, wormI, deadI, notwI, newpos ] = ClickAWorm( wboxi );
			if ( nclicks <= 0 )
				break;
			end
			wboxi.worms = [ wboxi.worms ; wboxi.deads(deadI, :) ; wboxi.notws(notwI, :) ; newpos ];
			wboxi.deads(deadI, :) = [];
			wboxi.notws(notwI, :) = [];
			wboxi.manual = true;
			MarkBox( wboxi, brdata );
		end
	elseif strcmpi(ans, 'E')
		disp('Pick Dead Worms');
		while 1
			[ nclicks, wormI, deadI, notwI, newpos ] = ClickAWorm( wboxi );
			if ( nclicks <= 0 )
				break;
			end
			wboxi.deads = [ wboxi.deads ; wboxi.worms(wormI, :) ; wboxi.notws(notwI, :) ; newpos ];
			wboxi.worms(wormI, :) = [];
			wboxi.notws(notwI, :) = [];
			wboxi.manual = true;
			MarkBox( wboxi, brdata );
		end
	elseif strcmpi(ans, 'X')
		disp('Pick Not Worms');
		while 1
			[ nclicks, wormI, deadI, notwI, newpos ] = ClickAWorm( wboxi );
			if ( nclicks <= 0 )
				break;
			end
			wboxi.notws = [ wboxi.notws ; wboxi.worms(wormI, :) ; wboxi.deads(deadI, :) ; newpos ];
			wboxi.worms(wormI, :) = [];
			wboxi.deads(deadI, :) = [];
			wboxi.manual = true;
			MarkBox( wboxi, brdata );
		end
	end
	cntdata.frames(wfno).wboxes(wwno) = wboxi;
end
return;

function [ xv, yv ] = SelectPoly();
xv = zeros(1, 0);
yv = zeros(1, 0);
hold on;
h = plot( [ NaN ], [ NaN ], '.-', 'Color', [ 0 0.5 0 ], 'LineWidth', 2 );
hold off;
while 1
	[ gx, gy ] = ginput(1);
	if isempty(gx)
		break;
	end
	xv(end+1) = gx(1);
	yv(end+1) = gy(1);
	set( h, 'XData', [ xv xv(1) ], 'YData', [ yv yv(1) ] );
end
delete(h);
return;
