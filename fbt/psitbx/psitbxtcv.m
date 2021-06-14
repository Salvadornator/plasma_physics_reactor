function y = psitbxtcv(shot,varargin)

%*PSITBXTCV	TCV LIUQE poloidal flux
% PSI = PSITBXTCV(SHOT[,T][,FORMAT][,SOURCE])
%	if SHOT is a string, interpreted as EQDSK file

tw = NaN; form = '*0';
if shot == -1 || shot >= 100000, from = 'ECON'; 
else                             from = 'LIUQE'; end
for k = 1:nargin-1
 switch class(varargin{k})
  case 'char'
   uppervar = upper(varargin{k});
   switch uppervar
    case {'*0' '-0' '+0' '01' 'FS' 'JPHI'}, form = uppervar;
    otherwise, from = varargin{k};
   end
  case 'double', tw = varargin{k};
 end
end

y = [];
switch upper(from)
 case {'FBTE','RAMP','FLAT','ECON'}, tree = 'PCS';     branch = '\MGAMS.PSITBX';
 case {'LIUQE','LIUQE2','LIUQE3'},   tree = 'RESULTS'; branch = '\PSITBX';
	%NB: tree 'results' is common to all LIUQEs. 'branch' is used only for LIUQE(1) in the 'FS' case.
 otherwise,                          tree = NaN;       disp('EFIT source assumed');
end

if isnan(tree),
 eqdskval=read_eqdsk(from,15,1);
% check COCOS
 [dummy,indr]=min(abs(eqdskval.rmesh-eqdskval.raxis));
 [dummy,indz]=min(abs(eqdskval.zmesh-eqdskval.zaxis));
% eqdskval.rjphi is a slight misnomer, it is in fact mu_0*R^2*p' + F*F'; this
%   will be off by 4*pi^2 when the COCOS index is in the wrong decade
 if eqdskval.rjphi(indr,indz)./eqdskval.gradshaf(indr,indz) > 2,
   eqdskval=read_eqdsk(from,2,1);
   eqdskval.psi=eqdskval.psi*2*pi;
   eqdskval.psimesh=eqdskval.psimesh*2*pi;
   eqdskval.psiaxis=eqdskval.psiaxis*2*pi;
   eqdskval.psiedge=eqdskval.psiedge*2*pi;
   eqdskval.pprime=eqdskval.pprime/(2*pi);
   eqdskval.FFprime=eqdskval.FFprime/(2*pi);
 end
 switch form
  case {'*0' '-0','+0','01'}
   x.data		= eqdskval.psi-eqdskval.psiedge;
   x.dim		= {eqdskval.rmesh',eqdskval.zmesh',0}';
   x.units	= 'T*m^2';
   if eqdskval.psiaxis-eqdskval.psiedge>=0,
     sip = +1;fip = '+0';
   else
     sip = -1;fip = '-0';
   end
  case 'FS'
   y=psitbxp2p(psitbxtcv(shot,tw,'*0',from),'FS');
   return
  case 'JPHI'
   psimeshabs = eqdskval.psiaxis+ ...
		eqdskval.psimesh.*(eqdskval.psiedge-eqdskval.psiaxis);
   pprime = interp1(psimeshabs,eqdskval.pprime,eqdskval.psi);
   FFprime = interp1(psimeshabs,eqdskval.FFprime,eqdskval.psi);
   rmesh = repmat(eqdskval.rmesh',1,size(pprime,2));
   x.data		= (pprime.*rmesh + FFprime./rmesh/(4e-7*pi))*2*pi;
   x.dim		= {eqdskval.rmesh',eqdskval.zmesh',0}';
   x.units	= 'A/m^2';
  otherwise, error('Invalid FORM')
 end   
else
 if isempty(mdsopen(tree,shot)), return, end
 if strcmp(tree,'PCS')
  ft = mdsdata('\MGAMS.DATA:AMAGIC') > 0;
  if strcmp(from,'ECON')
   if ft, from = 'RAMP'; else, from = 'FBTE'; end
  end
  if ~ft & ~strcmp(from,'FBTE')
   warning('RAMP or FLAT cannot be indentified in MGAMS')
   return
  end
 end

 ip = tdi('TCV_EQ("I_P",$1)',from);
 if isempty( ip.dim ), error( 'MDS data nodes are empty.' ); end
 t  = ip.dim{1};
 switch form
  case {'*0' '-0','+0','01'}
   x = tdi('TCV_EQ("PSI",$1)',from);
   sip = sign(mean(ip.data));
   if sip == 1, fip = '+0'; else, fip = '-0'; end
  case 'FS'
   if strcmp(from,'FBTE') & ft, tp = t([1:end/2 end/2:-1:1]); st = '$';
   else                         tp = t;                       st = '*';
   end
   % Calculate 'FS' from poloidal flux since precalculated values are not available in mds for LIUQE2 and LIUQE3
   if strcmp(from,'LIUQE2') | strcmp(from,'LIUQE3')
     y=psitbxp2p(psitbxtcv(shot,tw,'*0',from),'FS');
     return
   end
   x    =     tdi([branch ':AS[*,*,' st ']'],tp);
   rmag = mdsdata([branch ':RMAG  [' st ']'],tp);
   zmag = mdsdata([branch ':ZMAG  [' st ']'],tp);
   pmag = mdsdata([branch ':PSIMAG[' st ']'],tp);
   if st == '$', x.dim{3} = t; end
  case 'JPHI'
   x = tdi('TCV_EQ("J_TOR",$1)',from);
  otherwise, error('Invalid FORM')
 end

end

z = x.data;
t = x.dim{3};
x = x.dim(1:2);

if ~isnan(tw) & ~isempty(z)
 k = ifloor(round(t*1e6),round(tw(tw <= max(t))*1e6),1);
 k = k(k > 0); k(find(diff(k) == 0)) = [];
 if isempty( k ), error( [ 'Requested times are not in the data time interval: [', num2str( t(1) ), ', ', num2str( t(end) ), ']' ] ); end
 t = t(k);
 z = z(:,:,k);
 if strcmp(form,'FS')
  rmag = rmag(k); zmag = zmag(k); pmag = pmag(k);
 end
end

switch form
 case '*0'
  y = psitbxpsi(z,psitbxgrid('C','G',x),t,fip);
 case '+0'
  y = psitbxpsi(sip*z,psitbxgrid('C','G',x),t,'+0');
 case '-0'
  y = psitbxpsi(-sip*z,psitbxgrid('C','G',x),t,'-0');
 case '01'
  y = psitbxp2p(psitbxpsi(z,psitbxgrid('C','G',x),t,fip),'01');
 case 'FS'
  y = psitbxpsi(z,psitbxgrid('F','G',x),t,'FS',[rmag(:)';zmag(:)';pmag(:)']);
 case 'JPHI'
  y = psitbxfun(z,psitbxgrid('C','G',x),t);
end
