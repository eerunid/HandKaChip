function meta = ReadZeissMetaXML( xmlpath )
meta = struct;
meta.is_open = false;

xDoc = xmlread( xmlpath );

scaleX = NaN;
scaleY = NaN;
xScaling = xDoc.getElementsByTagName('Scaling');
if ( xScaling.getLength > 0 )
	xScaling = xScaling.item(xScaling.getLength-1);
	xDistances = xScaling.getElementsByTagName('Distance');
	for i = 1:xDistances.getLength
		xDISTi = xDistances.item(i-1);
		xATTRi = xDISTi.getAttributes;
		for j = 1:xATTRi.getLength
			attrname = xATTRi.item(j-1).getNodeName();
			attrvalue = char(xATTRi.item(j-1).getNodeValue());

			if ( strcmp( attrname, 'Id' ) )
				if ( strcmp( attrvalue, 'X' ) )
					scaleX = str2num( xDISTi.getElementsByTagName('Value').item(0).getTextContent() );
				end
				if ( strcmp( attrvalue, 'Y' ) )
					scaleY = str2num( xDISTi.getElementsByTagName('Value').item(0).getTextContent() );
				end
			end
		end
	end
end
if ( isnan(scaleX) || isnan(scaleY) )
	fprintf( 2, 'Could not retrieve Scaling Info (per Pixel).\n' );
end


xInfo = xDoc.getElementsByTagName('Information');
xInfo = xInfo.item(xInfo.getLength-1);
xImage = xInfo.getElementsByTagName('Image').item(0);
framew = str2num( xImage.getElementsByTagName('SizeX').item(0).getTextContent() );
frameh = str2num( xImage.getElementsByTagName('SizeY').item(0).getTextContent() );
nframes = 1;
xSizeT = xImage.getElementsByTagName('SizeT');
if ( xSizeT.getLength() > 0 )
	nframes = str2num( xSizeT.item(0).getTextContent() );
end
nchannels = str2num( xImage.getElementsByTagName('SizeC').item(0).getTextContent() );
nscenes = 1;
nstacks = 1;
xSizeS = xImage.getElementsByTagName('SizeS');
if ( xSizeS.getLength > 0 )
	nscenes = str2num( xSizeS.item(0).getTextContent() );
end
xSizeZ = xImage.getElementsByTagName('SizeZ');
if ( xSizeZ.getLength > 0 )
	nstacks = str2num( xSizeZ.item(0).getTextContent() );
end


xDimensions = xImage.getElementsByTagName('Dimensions').item(0);

% Channels
xChannels = xDimensions.getElementsByTagName('Channels').item(0).getElementsByTagName('Channel');
channels = [];
for i = 1:xChannels.getLength
	xCHi = xChannels.item(i-1);
	xATTRi = xCHi.getAttributes;
	for j = 1:xATTRi.getLength
		attrname = xATTRi.item(j-1).getNodeName();
		attrvalue = char(xATTRi.item(j-1).getNodeValue());

		if ( strcmp( attrname, 'Id' ) && strncmp( attrvalue, 'Channel:', 8 ) )
			chid = str2num( attrvalue(9:end) );
			channel = struct;
			channel.name = char(xCHi.getElementsByTagName('Fluor').item(0).getTextContent());
			binning  = char(xCHi.getElementsByTagName('DetectorSettings').item(0).getElementsByTagName('Binning').item(0).getTextContent());
			binning = sscanf( binning, '%d,%d' );
			channel.binningX = binning(1);
			channel.binningY = binning(2);
			channel.exposureT = str2num( xCHi.getElementsByTagName('ExposureTime').item(0).getTextContent() );
			channels = [ channels channel ];
		end
	end
end

% Time Series
xT = xDimensions.getElementsByTagName('T').item(0);
%tStart = = char( xT.getElementsByTagName('StartTime').item(0).getTextContent() );
%tStart = strrep( tStart, 'T', ' ' );
%tStart = strrep( tStart, 'Z', ' ' );
%jd = java.util.Date();
%tStart = datenum( tStart ) - jd.getTimezoneOffset()/(24*60);
%timestr = char( xT.getElementsByTagName('StartTime').item(0).getTextContent() );
timestr = char( xImage.getElementsByTagName('AcquisitionDateAndTime').item(0).getTextContent() );
timestr = strrep( timestr, 'T', ' ' );
timestr = strrep( timestr, 'Z', ' ' );
dotpos = strfind( timestr, '.' );
subsec = 0.0;
if ~isempty(dotpos)
	subsec = str2num(timestr(dotpos(end):end));
	timestr = timestr(1:(dotpos(1)-1));
end
tInc = 0;
xPositions = xT.getElementsByTagName('Positions');
if ( xPositions.getLength > 0 )
	xInterval = xPositions.item(0).getElementsByTagName('Interval');
	if ( xInterval.getLength > 0 )
		tInc = str2num( xInterval.item(0).getElementsByTagName('Increment').item(0).getTextContent() );
	end
end

% Z-Stacks
xZ = xDimensions.getElementsByTagName('Z');
stacks = [ 0 ];
if ( nstacks > 1 && xZ.getLength > 0 )
	xZ = xZ.item(0);
	xPositions = xZ.getElementsByTagName('Positions');
	if ( xPositions.getLength > 0 )
		xInterval = xPositions.item(0).getElementsByTagName('Interval');
		if ( xInterval.getLength > 0 )
			zStart = str2num( xInterval.item(0).getElementsByTagName('Start').item(0).getTextContent() );
			zInc = str2num( xInterval.item(0).getElementsByTagName('Increment').item(0).getTextContent() );
			stacks = zStart + zInc * ( 0:1:(nstacks-1) );
		end
	end
end


% Scenes
scenes = [];
posArrays = [];
xS = xDimensions.getElementsByTagName('S');
if ( xS.getLength > 0 )
	xS = xS.item(0);		
	xScenes = xS.getElementsByTagName('Scenes').item(0).getElementsByTagName('Scene');
	for i = 1:xScenes.getLength
		xScenesi = xScenes.item(i-1);
		xATTRi = xScenesi.getAttributes;

		scene = struct;
		for j = 1:xATTRi.getLength
			attrname = xATTRi.item(j-1).getNodeName();
			attrvalue =xATTRi.item(j-1).getNodeValue();

			scene = setfield( scene, char(attrname), char(attrvalue) );
		end

		xRegId = xScenesi.getElementsByTagName('RegionId');

		scene.Id = str2num(scene.Index);
		scene.RegId = str2num( [ 'uint64(' char(xRegId.item(0).getTextContent()) ')' ] );
		scene.PosName = scene.Name;
		scene.ArrayNo = 0;
		scenes = [ scenes scene ];
	end
	regIds = [ scenes.RegId ];

	xExperiment = xDoc.getElementsByTagName('Experiment').item(0);
	xSTRArrays_ = xExperiment.getElementsByTagName('SingleTileRegionArrays');

	totNumArrays = 0;
	for ii = 1:xSTRArrays_.getLength
		xSTRArrays = xSTRArrays_.item(ii-1).getElementsByTagName('SingleTileRegionArray');
		for i = 1:xSTRArrays.getLength
			xSTRArrayi = xSTRArrays.item(i-1);
			totNumArrays = totNumArrays + ( str2num( xSTRArrayi.getElementsByTagName('IsUsedForAcquisition').item(0).getTextContent() ) > 0 );
		end
	end
	
	if ( totNumArrays > 0 )
		for ii = 1:xSTRArrays_.getLength
			xSTRArrays = xSTRArrays_.item(ii-1).getElementsByTagName('SingleTileRegionArray');
		for i = 1:xSTRArrays.getLength
			xSTRArrayi = xSTRArrays.item(i-1);
			xATTRi = xSTRArrayi.getAttributes;

			if ( str2num( xSTRArrayi.getElementsByTagName('IsUsedForAcquisition').item(0).getTextContent() ) )

				array = struct;
				for j = 1:xATTRi.getLength
					attrname = xATTRi.item(j-1).getNodeName();
					attrvalue =xATTRi.item(j-1).getNodeValue();

					array = setfield( array, char(attrname), char(attrvalue) );
				end

				posArray = struct;
				posArray.name = array.Name;
				posArray.positions = cell(1, 0);
				posArray.sceneIds = [];

				xSTRs = xSTRArrayi.getElementsByTagName('SingleTileRegion');
				for ii = 1:xSTRs.getLength;
					xSTRii = xSTRs.item(ii-1);
					xATTRii = xSTRii.getAttributes;

					if ( str2num( xSTRii.getElementsByTagName('IsUsedForAcquisition').item(0).getTextContent() ) )
						posii = struct;
						for jj = 1:xATTRii.getLength
							attrname = xATTRii.item(jj-1).getNodeName();
							attrvalue =xATTRii.item(jj-1).getNodeValue();
	
							posii = setfield( posii, char(attrname), char(attrvalue) );
						end

						scenesii = find( str2num( [ 'uint64(' char(posii.Id) ')' ] ) == regIds );
						posii.scenes = reshape( scenesii, 1, [] );
						posArray.positions{end+1} = posii;
					end
				end
				if ~isempty(posArray.positions)
					posArrays = [ posArrays posArray ];
				end
			end
		end
		end
	else
		% No position arrays. Just single positions
		xSTRs_ = xExperiment.getElementsByTagName('SingleTileRegions');
		
		for ii = 1:xSTRs_.getLength
			if ~strcmp( xSTRs_.item(ii-1).getParentNode().getNodeName(), 'SingleTileRegionArray' )
				xSTRs = xSTRs_.item(ii-1).getElementsByTagName('SingleTileRegion');
				for i = 1:xSTRs.getLength
					xSTRi = xSTRs.item(i-1);
					xATTRi = xSTRi.getAttributes;

					if ( str2num( xSTRi.getElementsByTagName('IsUsedForAcquisition').item(0).getTextContent() ) )
	
						posi = struct;
						for j = 1:xATTRi.getLength
							attrname = xATTRi.item(j-1).getNodeName();
							attrvalue =xATTRi.item(j-1).getNodeValue();

							posi = setfield( posi, char(attrname), char(attrvalue) );
						end

						scenesi = find( str2num( [ 'uint64(' char(posi.Id) ')' ] ) == regIds );
						posi.scenes = reshape( scenesi, 1, [] );

						posArray = struct;
						posArray.name = posi.Name;
						posArray.positions = { posi };
						posArray.sceneIds = [];
						posArrays = [ posArrays posArray ];

					end
				end
			end
		end
	end
else
	% No Tile at all.
	scene = struct;
	scene.Index = '0';
	scene.Name = 'P0';
	scene.Id = 0;
	scene.PosName = 'P0';
	scene.ArrayNo = 0;
	scenes = [ scene ];

	pos = struct;
	pos.Name = 'P0';
	pos.scenes = 1;

	posArray = struct;
	posArray.name = 'P0';
	posArray.positions = { pos };
	posArray.sceneIds = [];
	posArrays = [ posArray ];
end

% Removed Incomplete Position Arrays
toremove = [];
scenes_available = 1:numel(scenes);
for i = 1:numel(posArrays)
	scenes_left = scenes_available;
	scenes_added = [];
	for k = 1:numel(posArrays(i).positions)
		scenesik = posArrays(i).positions{k}.scenes;
		while ~isempty(scenesik)
			indik = find( scenes_left == scenesik(1) );
			if ~isempty(indik)
				scenes_left(indik) = [];
				scenes_added(end+1) = scenesik(1);
				break;
			end
			scenesik(1) = [];
		end
	end
	if numel(scenes_added) >= numel(posArrays(i).positions)
		scenes_available = scenes_left;
		posArrays(i).sceneIds = reshape( [ scenes( scenes_added ).Id ], 1, [] );
	else
		toremove(end+1) = i;
	end
end
posArrays(toremove) = [];


meta = struct;
meta.is_open = true;
meta.px2umX = scaleX * 10^6;
meta.px2umY = scaleY * 10^6;
meta.frame2ms = tInc * 10^3;
meta.acqtime = datenum(timestr) + subsec/(24*3600) ;
meta.nframes = nframes;
meta.framew = framew;
meta.frameh = frameh;
meta.nchannels = nchannels;
meta.channels = channels;
meta.nstacks = nstacks;
meta.stacks = stacks;
meta.nscenes = nscenes;
meta.posArrays = posArrays;
meta.scenes = scenes;
return;
