function WriteZeissTData( path, fname, tdata );
fclose('all');

savpath = fullfile( path, [ fname '.tdata' ] );
fid = fopen( savpath, 'w' );
if fid == -1 
	fprintf( 2, 'TData[%s] open failed.\n', savpath );
	return;
end

fprintf( fid, 'Acquisition Time\r\n' );
fprintf( fid, '%.12f\r\n', tdata(1) );
fprintf( fid, 'Time (min)\r\n' );
for i = 1:numel(tdata)
	% LKSCMT: TData is in minutes as of 2014.12.23.
	fprintf( fid, '%4.0f\r\n', 24 * 60 * ( tdata(i) - tdata(1) ) );
end
fclose(fid);
return;
