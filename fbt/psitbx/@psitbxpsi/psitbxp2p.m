function psi = psitbxp2p(psi,mode,varargin)

% PSITBXPSI/PSITBXP2P

if strcmp( mode, '*0' )
    switch psi.format
        case {'-0' '+0'}
            mode = psi.format; %do nothing
        case {'01'}
            if ~isempty( psi.psimag )
                if mean( psi.psimag(:) ) < 0
                    mode = '-0';
                else
                    mode = '+0';
                end
            else
                error( 'Sign of flux cannot be deduced, .psimag is empty.' );
            end
        otherwise
            error( '''*0'' request is valid only for input objects of format ''01'', ''+0'' and ''-0''' )
    end
end

switch mode
 case {'-0' '+0'} 
  switch psi.format
   case {'-0' '+0'}
    if ~strcmp(mode,psi.format), 
        psi.psitbxfun = -psi.psitbxfun;
        psi.psimag = -psi.psimag; %JR 12.07.10
    end
   case '01' 
        psi.psitbxfun = (1 - psi.psitbxfun) .* ...
            repmat( abs( reshape( psi.psimag, 1, 1, [] ) ), [ size( psi.psitbxfun, [ 1, 2 ] ), 1] ) ...
            .* str2num( [ mode(1), '1' ] ); %NEW, JR 14.01.10
            %repmat(psi.psimag,[size(psi.psitbxfun.x(1)),size(psi.psitbxfun.x(2)),1]); %OLD, JR 13.01.10
        psi.psimag = abs( psi.psimag ) .* str2num( [ mode(1), '1' ] ); %JR 12.07.10
   case 'FS', error('"Flux-Surfaces" can not be converted')
  end
 case '01'
  psi = psitbxmag(psi);
 case 'FS'
  psi = psitbxfsd(psi,varargin{:});
 otherwise
  error('Invalid mode')
end
psi.format = mode;
