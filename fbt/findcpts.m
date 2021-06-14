function [rE,zE,fE,rS,zS,fS,rM,zM,fM,err]=findcpts(varargin)
% FINDCPTS - Find critical points (extrema and saddle points) in a 2D image.
%
%   The image has to be specified on an equidistant rectangular grid. The number
%   of grid location of the critical points is determined according to A. Kuijper, 
%   Pattern Recognition Letters, 25 (2004) 1665. The critical points are found
%   by Newton-Raphson iteration in the vicinity of the identified cell.
%
% CALL
%   >> [rE,zE,fE,rS,zS,fS,rM,zM,fM,err] = findcpts(psi[,ParameterName,ParameterValue]);
%
% INPUT
%   psi psitbxpsi or psitbxfun-object   Pol. flux on a clyindrical grid
% 
% PARAMETERS
%   itermax    {10}   Max. number of iterations in Newton-Raphson method to find critical
%                     points.
%   method     'Blom'
%   resolution        FLOAT: Cell sise.
%                     STRING: Starting with 'x' multiplier for the original grid size.
%                     The new cell size corresponds approx. to the min. distance between 
%                     critical points. 
%   tolerance  1e-6   Determine convergence if location changes by less than tolerance.
%   plot       [0,1]  Set to 1 to plot cells (default: 0).
%   verbose    [0,1]  Set to 1 to display additional information (default: 0).
%   debug      [0,1]  Set to 1 to debug (default:0).
%
% OUTPUT
%   rE  [max. # of extrema x #times]
%   zE
%   fE
%   rS  [max. # of saddle points x #times]
%   zS   
%   fS
%   rM  [max. # of monkey points x #times]  Refine grid to attempt a decomposition of 
%   zM                                      monkey points into two saddle points each.
%   fM
%   err        [0,1] 1 indicates that the algorithm likely missed a critical point. Best
%                    practice would be to re-run the case with a refined grid until err
%                    becomes zero.
%
% HISTORY
%   16-01-2013  Created by H. Reimerdes.
%   28-02-2013  Make sure that critical points do not leave the domain.                          
%   29-08-2013  Add possibility to refine grid (RESOLUTION), return monkey points and an 
%               error variable.
%   10-04-2015  Minor changes.
%

% See new_fincpts for a hexagonal grid.

  % --- Define default parameters ---
  p=struct( ...
    'itermax',    10, ...
    'method',     'Blom', ...
    'resolution',  [], ...        % Min. separation of critical points
    'tolerance',   1e-6, ...
    'plot',        0, ...
    'debug',       0, ...
    'verbose',     0 ...
    );
  
  switch class(varargin{1})
    case 'psitbxpsi'
      psi = varargin{1}.psitbxfun;
    case 'psitbxfun'
      psi = varargin{1};
    otherwise
      error('First input argument has to be a psitbxpsi or psitbxfun-object.')
  end
  
  % --- Overwrite default parameters ---
  for j=2:length(varargin)-1 
    if ischar(varargin{j}) && any(strcmpi(fieldnames(p),varargin{j}));
      p.(lower(varargin{j}))=varargin{j+1};
    end
  end

  if p.debug
    p.verbose=1;
    p.plot=1;
  end
 
  if p.verbose, fprintf('Starting FINDCPTS ...\n'); end
  
  err = 0;
  
  % Dimensions of original flux map (Nr,Nz,Nt)
  [Nr,Nz,~] = size(psi);
  if isempty(psi.t), Nt = 1; else Nt=length(psi.t); end
   
  % Grid 
  g = psi.grid;
  
  % Grid spacing (has to be equidistant!)
  if ~strcmp(g.type(1),'C') || ~strcmp(g.storage(1),'G')
    error('Grid points have to be stored on a grid in cylindrical coordinates');
  end

  r = psi.grid.x{1}; rmin=min(r); rmax=max(r);
  z = psi.grid.x{2}; zmin=min(z); zmax=max(z);
  Dr = mean(diff(g.x{1}));
  Dz = mean(diff(g.x{2}));
  
  if p.verbose, fprintf('\tGrid resolution: Delta R=%0.3f, Delta Z=%0.3f\n',Dr,Dz); end 
  
  %
  switch lower(p.method)
    case 'blom'
     
      if ~isempty(p.resolution) && ( ischar(p.resolution) || Dr > p.resolution || Dz > p.resolution )
        
        if ischar(p.resolution)
          Nr0=round((Nr-1)*sscanf(p.resolution,'x%g'))+1;
          Nz0=round((Nz-1)*sscanf(p.resolution,'x%g'))+1;
        else
          Nr0=round((rmax-rmin)/p.resolution)+1;
          Nz0=round((zmax-zmin)/p.resolution)+1;
        end
        
        g0 = psitbxgrid('C','G',{linspace(rmin,rmax,Nr0), ...
                                 linspace(zmin,zmax,Nz0),NaN});
        Dr0 = (rmax-rmin)/(Nr0-1); 
        Dz0 = (zmax-zmin)/(Nz0-1);  
        
        if p.verbose
          fprintf('\tRefine grid from %0d to %0d radial and %0d to %0d vertical grid points.\n', ...
            Nr,Nr0,Nz,Nz0);
        end
      else
        g0 = g;
        Dr0 = Dr;
        Dz0 = Dz;   
        Nr0 = Nr;
        Nz0 = Nz;
      end
    
  end
  
  % Initialise variables for extrema, saddle and monkey points
  rE = []; zE = []; fE = [];
  rS = []; zS = []; fS = [];
  rM = []; zM = []; fM = [];
  
  for j=1:Nt % Loop over times 
    
    if isempty(psi.t), t=[]; else t = psi.t(j); end
    
    f = psitbxfun(psi.x(:,:,j),g,t); 
    
    % --- Investigate neighbourhoods of each pixel by number of sign changes ---
    nbh = nan(Nr0,Nz0);
    
    % Function values on the finer grid 
    f0 = psitbxf2f(f,g0);
    
    % Indices of neighbourhood
    dk = [ 0,-1,-1, 0,+1,+1, 0];  % 1st dim, i.e. row index - clockwise, closed
  
    for k=2:Nr0-1 % Loop over rows
      
      % 2nd dim, i.e. column index, depends on row number (Blom translation)
      dl=[-1, -mod(k,2), -mod(k,2)+1, 1,  -mod(k,2)+1, -mod(k,2), -1];
 
      for l=2:Nz0-1
        dfkl = f0.x(sub2ind([Nr0,Nz0],k+dk,l+dl)) > f0.x(k,l);  % 1 for + and 0 for -
        nbh(k,l)=length(find(diff(dfkl)));
      end 
      
    end % of loop over rows
  
    if p.verbose
      fprintf('\tAnalysing time slice %0d/%0d:\n',j,Nt);  
      fprintf('\tZero sign changes: %d\n',length(find(nbh==0)));
      fprintf('\tTwo sign changes: %d\n',length(find(nbh==2)));
      fprintf('\tFour sign changes: %d\n',length(find(nbh==4)));
      fprintf('\tSix sign changes: %d\n',length(find(nbh==6))); 
    end
    
    % ================================ Critical points ===================================
    
   
    % Critical points for j-th time step   
    rC0=[]; zC0=[]; tC0=[];                         % Cell
    rCj = []; zCj = []; fCj = []; cCj = []; NCj=0;  % Iterated position

    [ir,iz]=find(nbh==0 | nbh==4 | nbh==6);
    if ~isempty(ir)
      rC0 = g0.x{1}(ir);
      zC0 = g0.x{2}(iz);
      tC0 = nbh(sub2ind([Nr0,Nz0],ir,iz))';
     
      % Set up variables for iteration
      NCj = length(rC0);    
      [rCj,zCj,fCj] = deal(nan(1,NCj));
      cCj = zeros(1,NCj);

      % Calculate derivatives
      fr = psitbxf2f(f,[1 0]); % First derivatives
      fz = psitbxf2f(f,[0 1]);      
      frr = psitbxf2f(fr,[1 0]); % Second derivatives
      frz = psitbxf2f(fr,[0 1]);
      fzz = psitbxf2f(fz,[0 1]);
    end
    
    if p.verbose, fprintf('\tLocating %d critical points\n',length(ir)); end
    
    
    for k=1:NCj  % Loop over critical points
      
      switch tC0(k)
        case 0
          tC0k='extremum';
        case 4
          tC0k='saddle point';
        case 6
          tC0k='monkey point';
      end
        
      if p.verbose
        fprintf('\t%0d: Refine location of %s in the vicinity of (R0,Z0)=(%0.3f,%0.3f)m\n', ...
          k,tC0k,rC0(k),zC0(k));
      end
      gCk = psitbxgrid('C','Points',{rC0(k),zC0(k),NaN});
         
      converged=0;  iter=0;
       
      while ~converged && iter<p.itermax    % Newton-Raphson iteration

        f0 = psitbxf2f(f,gCk);
        fr0 = psitbxf2f(fr,gCk);
        fz0 = psitbxf2f(fz,gCk);       
        frr0 = psitbxf2f(frr,gCk);
        frz0 = psitbxf2f(frz,gCk);
        fzz0 = psitbxf2f(fzz,gCk);
            
        dr = (fz0.*frz0 - fr0.*fzz0)./(frr0.*fzz0 - frz0.^2);
        dz = (fr0.*frz0 - fz0.*frr0)./(frr0.*fzz0 - frz0.^2);
        
        % New estimate of critical points
        gCk_new=psitbxgrid('C','Points',{gCk.x{1} + dr.x, gCk.x{2} + dz.x, NaN});
        
        if gCk_new.x{1}<f.grid.x{1}(1) || gCk_new.x{1}>f.grid.x{1}(end) ...
         || gCk_new.x{2}<f.grid.x{2}(1) || gCk_new.x{2}>f.grid.x{2}(end)     
          fprintf('\t\tWARNING: %s with (R0,Z0)=(%0.3f,%0.3f)m leaves computational grid.\n', ...
            tC0k,rC0(k),zC0(k));
          break
        else
           % Prepare next iteration
          iter=iter+1;
          gCk=gCk_new;
          
          % Check whether extremum leaves cell
          if abs(gCk.x{1}-rC0(k))>Dr || abs(gCk.x{2}-zC0(k))>Dz
            fprintf('\t\tWARNING: %s at (R,Z)=(%0.3f,%0.3f)m leaves original cell (%0.3f,%0.3f)m -> consider decreasing the resolution.\n', ...
              tC0k,gCk.x{1},gCk.x{2},rC0(k),zC0(k));
             % Debug
            if p.debug
              figure
              if ~exist('gh','var')
                gh=psitbxgrid('C','G',{linspace(f.grid.x{1}(1),f.grid.x{1}(end),501),linspace(f.grid.x{2}(1),f.grid.x{2}(end),1001),NaN});
                fh=psitbxf2f(f,gh);
              end
              contour(fh.grid.x{1},fh.grid.x{2},fh.x',201)
              hold on
              plot(g0,'+r')
              plot(f.grid,'ok')
              plot([rC0(k),gCk.x{1}],[zC0(k),gCk.x{2}],'-or')
              set(gca,'xlim',rC0(k)+[-0.1 0.1],'ylim',zC0(k)+[-0.1,0.1])
              keyboard; 
            end
          end
          
          % Check step size
          if sqrt(dr.x.^2+dz.x.^2) < p.tolerance
            converged = 1;
            if p.verbose, fprintf('\t\t%s location converged after %0d iterations -> (%0.3f,%0.3f)m\n',tC0k,iter,gCk.x{1},gCk.x{2}); end
          else
            if (iter == p.itermax && p.verbose)
              fprintf('\t\tNewton-Raphson method exceeds max. number of iterations.\n')
            end
          end 
        end  
      end % of iteration (while loop)
      
      fCk=psitbxf2f(f,gCk);
      
      rCj(k) = gCk.x{1};
      zCj(k) = gCk.x{2};
      fCj(k) = fCk.x;
      cCj(k) = converged;
      
    end % of loop over critical points (k)  
    
    for k=1:NCj-1
      ik = find(cCj([k+1:NCj]) & abs(rCj([k+1:NCj])-rCj(k))<Dr0/2 & abs(zCj([k+1:NCj])-zCj(k))<Dz0/2);
      if cCj(k) && ~isempty(ik)
        ik=[k, k+ik];
        fprintf('\t\tWARNING: Found identical points: %s -> refine grid to locate all critical points.\n',num2str(ik));
        ii = imin(sqrt((rCj(ik)-rC0(ik)).^2+(zCj(ik)-zC0(ik)).^2));
        cCj(ik)=0;
        cCj(ik(ii))=1; % Keep the point that is closest to its original cell centre
        err = 1;    % Signal that a critical point is missed. 
      end
    end
    
    % --- Assign to variables according to type ---
    iE = find(tC0==0 & cCj);
    iS = find(tC0==4 & cCj);
    iM = find(tC0==6 & cCj);
        
    fEj = fCj(iE);  rEj = rCj(iE); zEj = zCj(iE);
    fSj = fCj(iS);  rSj = rCj(iS); zSj = zCj(iS);
    fMj = fCj(iM);  rMj = rCj(iM); zMj = zCj(iM);   
    
    % --- Update output variables ---
    rE = addcol(rE,rEj);  zE = addcol(zE,zEj);  fE = addcol(fE,fEj);
    rS = addcol(rS,rSj);  zS = addcol(zS,zSj);  fS = addcol(fS,fSj);
    rM = addcol(rM,rMj);  zM = addcol(zM,zMj);  fM = addcol(fM,fMj);
    
    % -- Plot result ---
    if p.plot
      if ~exist('hav','var') || ~ishandle(hav)
        figure('Name','FINDCPTS','Position',[0 0 380,640]);
        hav=subplot(1,1,1);
      else
        axes(hav); 
        cla; 
      end
      rm = 0.5*(g0.x{1}(1:end-1)+g0.x{1}(2:end));
      zm = 0.5*(g0.x{2}(1:end-1)+g0.x{2}(2:end)); 
      pcolor(hav,rm,zm,nbh(2:end,2:end)')
      title({sprintf('Number of sign changes for t=%0.3fs',psi.t(j)), ...
          '0: extremum, 2: normal, 4: saddle point, 6: monkey point'})
      set(hav,'clim',[0 6])
      axis equal
      colorbar
      hold on
      
      % Re-calculate derivatives in case it hadn't been done before
      fr = psitbxf2f(f,[1 0]);
      fz = psitbxf2f(f,[0 1]);
    
      contour(hav,g.x{1},g.x{2},fr.x',[0 0],'c')
      contour(hav,g.x{1},g.x{2},fz.x',[0 0],'m')
	 
      plot(hav,rEj,zEj,'w+')
      plot(hav,rSj,zSj,'wx')
      fprintf('\tType ''return'' to continue.\n')
      keyboard
    end

  end % of loop over times
  
  if p.verbose, fprintf('Exciting FINDCPTS normally.\n'); end
   
end % of function findcpts


function out=addcol(array,col)
    
  dim=size(array);
       
  if length(col)<=dim(1)
    out = [array ,[col(:); nan(dim(1)-length(col),1)]];
  else
    out = [[array; nan(length(col)-dim(1),dim(2))],col(:)];
  end
    
end     