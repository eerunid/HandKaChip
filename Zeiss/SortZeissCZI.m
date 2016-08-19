function czi = SortZeissCZI( czipath )
czi = struct;
czi.is_open = false;

fid = fopen( czipath, 'r' );
if ( fid == -1 )
	fprintf( 2, 'Failed to open CZI file[%s].\n', czipath );
	return;
end

tdata = cell( 1, 0 );
opword = '<METADATA>';
clword = '</METADATA>';
oplen = length(opword);
cllen = length(clword);

global opword0 clword0 oplen0 cllen0
opword0 = '<AcquisitionTime>';
clword0 = '</AcquisitionTime>';
oplen0 = length(opword0);
cllen0 = length(clword0);

data = [];
segsize = 10000;
while ~feof(fid)
	temp = fread( fid, segsize, '*char' )';
	data = [ data temp ];

	is_open = false;
	sind = strfind( data, opword );
	eind = [ sind(2:end)-1 length(data) ];
	nfound = numel(sind);
	for i = 1:nfound
		datai = data(sind(i):eind(i));
		indi = strfind( datai, clword );
		if numel(indi) > 1
			fprintf( 2, 'Too many closes: [%d] closes\n', numel(indi) );
		end

		if isempty(indi)
			if i >= nfound
				is_open = true;
			else
				fprintf( 2, 'Unmatched parenthesis\n');
			end
		else
			tdatai = ProcessMetadata( datai(1:(indi(1)+cllen-1)) );
			if ~isempty(tdatai)
				tdata{end+1} = tdatai;
			end
		end
	end

	if is_open
		data = data(sind(end):end);
	elseif length(data) > oplen
		data = data((end-oplen+1):end);
	else
		data = [];
	end
end
fclose(fid);

czi.tdata_datenum = zeros( size(tdata) );
czi.tdata_subsec = zeros( size(tdata) );
for i = 1:numel(tdata)
	timestr = tdata{i};
	timestr = strrep( timestr, 'T', ' ' );
	timestr = strrep( timestr, 'Z', ' ' );
	dotpos = strfind( timestr, '.' );
	subsec = 0.0;
	if ~isempty(dotpos)
		subsec = str2num(timestr(dotpos(end):end));
		timestr = timestr(1:(dotpos(1)-1));
	end
	czi.tdata_datenum(i) = datenum(timestr);
	czi.tdata_subsec(i) = subsec;
end
czi.is_open = true;
return;

function timestamp = ProcessMetadata( data );
timestamp = [];
global opword0 clword0 oplen0 cllen0

sind = strfind( data, opword0 );
eind = [ sind(2:end)-1 length(data) ];
nfound = numel(sind);
for i = 1:nfound
	datai = data((sind(i)+oplen0):eind(i));
	indi = strfind( datai, clword0 );
	if ~isempty(indi)
		timestamp = datai(1:(indi(1)-1));
		return;
	end
end
return;
