function status = WriteZeissHeader( path, fname, acqtime, px2um, frame2ms, nframes, channels, stacks, patches );
status = -1;

hdrpath = fullfile( path , [ fname '.dat' ] );
fid = fopen( hdrpath, 'w' );
if ( fid == -1 )
	fprintf( 2, 'File open failed: %s \n', hdrpath );
	return;
end

fprintf( fid, 'Filename\r\n' );
fprintf( fid, '%s\r\n', fname );

fprintf( fid, 'Acquisition Time\r\n');
fprintf( fid, '%.12f\r\n', acqtime );

fprintf( fid, 'Length Per Pixel (um/px)\r\n');
fprintf( fid, '%.12g\r\n', px2um );

fprintf( fid, 'Frame Rate (ms/frame)\r\n');
fprintf( fid, '%.12g\r\n', frame2ms );

fprintf( fid, 'Number of Frames\r\n');
fprintf( fid, '%d\r\n', nframes );

fprintf( fid, 'Number of Channels\r\n');
fprintf( fid, '%d\r\n', numel(channels) );
for i = 1:numel(channels)
	fprintf( fid, '"%s" %d %d %d\r\n', ...
				channels(i).name, channels(i).binX, channels(i).binY, channels(i).expT );
end

fprintf( fid, 'Number of Z-Stacks\r\n');
fprintf( fid, '%d\r\n', numel(stacks) );
for i = 1:numel(stacks)
	fprintf( fid, '"%s" %.3f\r\n', ...
				stacks(i).name, stacks(i).zpos );
end

fprintf( fid, 'Number of Patches\r\n');
fprintf( fid, '%d\r\n', numel(patches) );
for i = 1:numel(patches)
	fprintf( fid, '"%s" %d %d %d %d\r\n', ...
				patches(i).name, round( [ patches(i).l patches(i).r patches(i).b patches(i).t ] ) );
end

fclose(fid);
status = 0;
return;
