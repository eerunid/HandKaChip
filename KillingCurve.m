function KillingCurve;
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 15);
clear all;
fclose('all');
workpath = input('Directory : ', 's');
[ tmpp, workname ] = fileparts( workpath );

mon = get(0, 'MonitorPositions' );
mon(:, [ 2 4 ] ) = -mon( :, [ 4 2 ] ) + max( mon(:, 4) ) + 1;
mx = mon(end, 1);
my = mon(end, 2);
mw = mon(end, 3)-mon(end, 1)+1;
mh = mon(end, 4)-mon(end, 2)+1;

figure;
clf;
%figx = mx + ( mod( gcf-1, 4 ) ) * mw/4;
%figy = my + mh - ceil( (mod( gcf-1, 8 )+1)/4 ) * mh/2;
%set( gcf, 'units', 'pixels', 'outerposition', [ figx figy mw/4 mh/2 ] );

global data
data = [];
FindAllKACntIn(workpath);
ndata = numel(data);

workpathlen = length(workpath)+1;
dataname = cell( 1, 0 );
for i = 1:ndata
	[ pathi, fnamei ] = fileparts( data(i).path );
	KAi = LoadKACnt( data(i).path );
	fnamei = fnamei(1:end-length('.KACnt'));
	tdatai = ReadTData( pathi, fnamei );
	%tdatai = 24 * ( tdatai - tdatai(KAi.selected(1)) );
	%tdatai = 24 * ( tdatai - tdatai(1) );
	if isempty(tdatai)
		tdatai = KAi.selected;
	end

	KAi.nframes = numel(KAi.selected);
	KAi.nworms = KAi.nworms(KAi.selected);
	KAi.tdata = tdatai(KAi.selected);

	%data(i).name = pathi(spos:end);
	spos = [ strfind( fnamei, ' ' ) 0 ] + 1;
	data(i).name = fnamei(spos:end);
	data(i).KA = KAi;
	dataname{i} = data(i).name;
end
[ dataname, sortI ] = sort( dataname );
data = data(sortI);
[ udataname, uniqueI ] = unique(dataname);

errtol = 0;
%hsvn = hsv(ndata);
hsvn = hsv(numel(uniqueI));
dataname = udataname;
cmapid = 1:ndata;
cnameid = 1:ndata;
for i = 1:numel(uniqueI)
	u = uniqueI(i);
	cmapid(u:end) = i;
	cnameid(u:end) = i;
end
cmap = zeros( ndata, 3 );
cname = cell( 1, ndata );
for i = 1:ndata
	%cmap(i, :) = hsvn( mod(i-1, 4)+1, : );
	%cname{i} = dataname{ mod(i-1, 4)+1 };
	cmap(i, :) = hsvn( cmapid(i), : );
	cname{i} = dataname{ cnameid(i) };
end

for i = 1:ndata
	ntotal = max(data(i).KA.nworms);
	data(i).KA.ntotal = ntotal;
end

hold on;
for u = reshape( uniqueI, 1, [] )
	plot( 0, 0, '-', 'Color', cmap(u, :), 'LineWidth', 2, 'DisplayName', cname{u} );
end
legend( 'show' );
for i = 1:ndata
	plot( data(i).KA.tdata, data(i).KA.nworms/data(i).KA.ntotal * 100, ...
			'o-', 'Color', cmap(i, :), ...
			'LineWidth', 2, 'DisplayName', cname{i} );
end
hCurrent = plot( [ NaN ], [ NaN ], 'o-', 'Color', [ 0 0 0 ], 'LineWidth', 2 );
hold off;
xlabel( 'Time (h)' );
ylabel( '% Survived' );
title( workname );

i = 0;
while 1
	i = max( 0, i );
	i = min( i, ndata );

	if i > 0
		fprintf( 1, '%s \n', data(i).path(workpathlen:end) );
		set( hCurrent, 'XData', data(i).KA.tdata, ...
						'YData', data(i).KA.nworms/data(i).KA.ntotal * 100 );
	else
		set( hCurrent, 'XData', [ NaN ], 'YData', [ NaN ] );
	end

	ans = input('[Q]uit ? ', 's');
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
		temp = round( input( sprintf( 'Data No. ? [ 0 - %d ] ', ndata ) ) );
		if ( ~isempty(temp) )
			i = temp;
		end
	end
end

return;

function FindAllKACntIn(path);
global data

if isempty(path)
	return;
end
files = dir(path);
for i = 1:length(files)
	if ( files(i).isdir )
		if ( ~strcmp(files(i).name, '.') && ~strcmp(files(i).name, '..') )
			FindAllKACntIn( fullfile( path, files(i).name ) );
		end
	end
end
files = dir( fullfile( path, '*.KACnt.mat' ) );
for i = 1:length(files)
	datum = struct;
	datum.path = fullfile( path, files(i).name );
	datum.name = '';
	datum.KA = struct;
	data = [ data datum ];
end
return;

function KA = LoadKACnt( path );
% Load .KACnt.mat
load( path, '-mat', 'cntdata' );

nframes = numel(cntdata.frames);
nworms = zeros( 1, nframes );
for i = 1:nframes
	wboxes = cntdata.frames(i).wboxes;
	wormsi = cat( 1, wboxes.worms );
	nworms(i) = size( wormsi, 1 );
end

KA = struct;
KA.selected = cntdata.selected;
KA.nframes = nframes;
KA.nworms = nworms;
return;

function tdata = ReadTData( path, fname );
tdata = [];
if exist( fullfile( path, [ fname '.tdata' ] ), 'file' ) > 0
	data = importdata( fullfile( path, [ fname '.tdata' ] ), '', 1 );
	acqtime = data.data;
	data = importdata( fullfile( path, [ fname '.tdata' ] ), '', 3 );
	tdata = data.data ;
	% LKSCMT: TData is in minutes as of 2014.12.23.
	tdata = reshape( tdata / 60, 1, [] );	
end
return;
