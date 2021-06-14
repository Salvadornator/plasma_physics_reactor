function pobj = eqdsk2psitbx(varargin)
% EQDSK2PSITBX
%   Transform an eqdsk equilibrium into a PsiTbx poloidal flux object.
%
% CALL
%   >> psi=eqdsk2psitbx(eqdsk_filename);
% or
%   >> psi=eqdsk2psitbx(eq_struct);
%
% INPUT
%   eqdsk_filename CHAR 
%   eq_struct      STRUCT   Output of read_eqdsk.m
%
% OUTPUT
%   psi   PsiTbx poloidal flux object
%
% HISTORY
%   11-08-2014  HR  New
%

  switch class(varargin{1})
    case 'char'
      eq = read_eqdsk_fast(varargin{1}); % Returns flux per radian
    case 'struct'
      eq = varargin{1};
  end

  gridp = psitbxgrid('Cylindrical','Grid',{eq.rmesh(:)';eq.zmesh(:)';NaN});
  x = eq.psi - eq.psiedge - 0.005*(eq.psiaxis-eq.psiedge);  % Make sure that psi=0 contour is closed
  if (eq.psiaxis>eq.psiedge), form='+0'; else form='-0'; end 
  t = -1.0;
  pobj = psitbxpsi(x.*(2*pi),gridp,t,form);

end
