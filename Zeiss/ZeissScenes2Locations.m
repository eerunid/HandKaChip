function [ nscenes, scenes, locations ] = ZeissScenes2Locations( posArrays, scenes );
nscenes = numel(scenes);

locations = [];
if ( isempty(posArrays) )
	sind = 1:numel(scenes);

	location = struct;
	location.name = '';
	location.patches = [];
	locations = repmat( location, 1, numel(sind) );
	ndigit = 1 + floor( log10( numel(sind) ) ) ;
	for lcno = 1:numel(sind)
		lcstr = num2str(lcno);
		locations(lcno).name = [ 'P' repmat( '0', 1, ndigit-length(lcstr) ) lcstr ];

		patches = repmat( struct, 1, 1 );
		patches(1).name = '';
		% In MATLAB, numbering starts from 1, not 0
		patches(1).l = scenes(sind(lcno)).StartX + 1;
		patches(1).r = scenes(sind(lcno)).StartX + scenes(sind(lcno)).SizeX -1 + 1;
		patches(1).b = scenes(sind(lcno)).StartY + 1;
		patches(1).t = scenes(sind(lcno)).StartY + scenes(sind(lcno)).SizeY -1 + 1;
		locations(lcno).patches = patches;

		scenes(sind(lcno)).lcno = lcno;
		scenes(sind(lcno)).pcno = 1;
	end
else
	for ano = 1:numel(posArrays)
		location = struct;
		location.name = posArrays(ano).name;
		location.patches = [];
		locations = [ locations location ];
		lcno = numel(locations);

		[ C, ia, ib ] = intersect( [ scenes.StartS ] , posArrays(ano).sceneIds );
		sind = ia;

		if ( isempty(sind) )
			locations(end) = [];
		else
			ndigit = 1 + floor( log10( numel(sind) ) ) ;
			patches = repmat( struct, 1, numel(sind) );
			for i = 1:numel(sind)
				pcstr = num2str(i);
				patches(i).name = [ repmat( '0', 1, ndigit-length(pcstr) ) pcstr ];
				% In MATLAB, numbering starts from 1, not 0
				patches(i).l = scenes(sind(i)).StartX + 1;
				patches(i).r = scenes(sind(i)).StartX + scenes(sind(i)).SizeX -1 + 1;
				patches(i).b = scenes(sind(i)).StartY + 1;
				patches(i).t = scenes(sind(i)).StartY + scenes(sind(i)).SizeY -1 + 1;
				scenes(sind(i)).lcno = lcno;
				scenes(sind(i)).pcno = i;
			end
			locations(lcno).patches = patches;
			if ( numel(patches) <= 1 )
				locations(lcno).patches(1).name = '';
			end
		end
	end
%	if ( numel(locations) <= 1 )
%		locations(1).name = '';
%	end
end

return;
