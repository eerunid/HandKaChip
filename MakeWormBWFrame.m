function BW = MakeWormBWFrame( frame, param );
sz = size(frame);
framew = sz(2);
frameh = sz(1);

%	[ X, Y ] = meshgrid( -param.wormR:param.wormR, -param.wormR:param.wormR );
%	mask = double( X.^2 + Y.^2 <= param.wormR^2 );
%	BWmask.mask = mask / sum(mask(:));
intR = ceil( param.wormR );
extR = ceil( param.bgndR );
[ X, Y ] = meshgrid( -extR:extR, -extR:extR );
intI = find( X.^2 + Y.^2 <= param.wormR^2 );
extI = setdiff( find( X.^2 + Y.^2 <= param.bgndR^2 ), intI );
mask = zeros( size(X+Y) );
mask(extI) = -1/numel(extI);
mask(intI) = +1/numel(intI);

cframe = conv2( frame, mask, 'same' );
%sdata = sort(frame(:), 'ascend');
%diff = sdata(ceil(0.8*end))-sdata(ceil(0.1*end));
BW = cframe < -1.0*std(frame(:));



% Expand BW
BW = double(BW);
exBW = BW;

[ X, Y ] = meshgrid( -intR:intR, -intR:intR );
intI = find( X.^2 + Y.^2 <= param.wormR^2 );
mask = zeros( size(X+Y) );
mask(intI) = +1/numel(intI);
cframe = conv2( frame, mask, 'same' );

CC = bwconncomp(BW);
S = regionprops(CC, 'Perimeter');
for i = 1:CC.NumObjects
	if ( S(i).Perimeter > param.min_wormP )
		indi = CC.PixelIdxList{i};
		sdatai = sort( cframe(indi), 'ascend' );
		cutoffi = sdatai( max( 1, ceil( 0.95*end ) ) );

		iindi = indi;
		while ~isempty(iindi)
			wi = ceil(iindi/frameh);
			hi = iindi - frameh*(wi-1);

			xindi = [];
			xindi = union( xindi, iindi(find(hi>1))-1 );
			xindi = union( xindi, iindi(find(hi<frameh))+1 );
			xindi = union( xindi, iindi(find(wi>1))-frameh );
			xindi = union( xindi, iindi(find(wi<framew))+frameh );
			xindi( find( exBW(xindi) > 0 ) ) = [];
			xindi( find( cframe(xindi) > cutoffi ) ) = [];

			if isempty(xindi)
				% No more pixel to extend
				break;
			end
			indi = union( indi, xindi );
			iindi = xindi;

			if ( numel(indi) > 1.5 * sum( BW(indi) ) )
				% Too many extended pixels
				break;
			end
			exBW(xindi) = max( 0.5, exBW(xindi) );
		end
	end
end
BW = exBW;
return;
