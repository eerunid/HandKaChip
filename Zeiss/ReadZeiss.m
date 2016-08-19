function [ frame, files ] = ReadZeiss( zeiss, chno, stno, frameno );
frame = [];
files = cell( 1, 0 );
if ( zeiss.is_open )
	if ( zeiss.bmask(frameno) )
		frame = zeiss.bframes{frameno}( :, :, chno, stno );
	else
		if ( ~isempty(zeiss.frames) )
			framefname = zeiss.frames{frameno};
		else
			ndigit = 1 + floor( log10( zeiss.nframes ) );
			framestr = num2str(frameno);
			framestr = [ repmat( '0', 1, ndigit-length(framestr) ) framestr ];
			framefname = [ zeiss.fname ' t' framestr '.tif' ];
		end

		if ( nargout <= 1 )
			frame = zeros( zeiss.frameh, zeiss.framew );
		end
		files = cell( 1, zeiss.npatches );
		% LKSCMT : Why descend?
		% The earlier patch is taken earlier, which means less photobleaching.
		for i = zeiss.npatches:-1:1
			patch = zeiss.patches(i);
			files{i} = fullfile( zeiss.path, zeiss.channels(chno).name, zeiss.stacks(stno).name, patch.name, framefname );

			if ( nargout <= 1 )
				temp = dir( files{i} );
				%if ( temp.bytes > 0 )
				if ( ~isempty(temp) && temp.bytes > 0 )
					pframe = double( imread( files{i} ) );
					frame( patch.b:patch.t, patch.l:patch.r ) = pframe;
				end
			end
		end
	end
end
return;
