function info = ReadZeissInfoXML( xmlpath )
info = struct;
info.is_open = false;

xDoc = xmlread( xmlpath );
xBounds = xDoc.getElementsByTagName('Bounds');

scenes = [];
for i = 1:xBounds.getLength
	xBoundsi = xBounds.item(i-1);
	xATTRi = xBoundsi.getAttributes;

	scene = struct;
	scene.StartZ = 0;
	scene.StartC = 0;
	scene.StartT = 0;
	scene.StartS = 0;
	for j = 1:xATTRi.getLength()
		attrname = xATTRi.item(j-1).getNodeName();
		attrvalue =xATTRi.item(j-1).getNodeValue();

		scene = setfield( scene, char(attrname), str2num(attrvalue) );
	end

%	if ( isempty(scenes) )
%		scenes = [ scene ];
%	elseif ~ismember( scene.StartS, [ scenes(:).StartS ] )
%		scenes = [ scenes scene ];
%	else
%		break;
%	end
	scenes = [ scenes scene ];
end

info.is_open = true;
info.nscenes = numel(scenes);
info.scenes = scenes;
return;
