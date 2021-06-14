function params = fbt(params)

  % --- Start running FBT ---
  fprintf('\n   --- Starting FBT ...\n'); 
  
  % --- Get FBT parameters structure ---
  %params = get_fbt_params;

  % --- Checking if Nb = length(Rb) = length(Zb) ---
  if any(params.Nb ~= max(size([params.Rb; params.Zb])))
      fprintf('   --- Warning: length of Nb, Rb and Zb must match\n')
      fprintf('   --- Stopping FBT\n')
      return
  end
  
  % --- Loading kinetic profiles ---
  params = getprofs(params);

  % --- Initial plasma current density distribution ---
  params = j0_guess(params,0);
  
  err = 1;
  imsg = 0;
  while err > params.err && params.iter < params.Niter
      
      params.iter = params.iter + 1;
      fprintf('   --- Iteration %2.0u\n',params.iter)
      
      % --- Updating plasma current density ---
      if params.iter > 1
        params = new_jplasma(params,0);
      end
      
      Ia = params.Ia;
      % --- Calculating plasma poloidal flux distribution ---
      params = j2psi(params,0);

      % --- Calculating Ohmic and vessel currents ---
      if params.iter > 1
          params = get_iv_ioh(params);
      end
      
      % --- Interpolating poloidal flux at specific points ---
      params = get_psivalues(params);
      
      % --- Calculating minimum active coil currents ---
      params = calc_currents(params);

      % --- Solving Grad-Shafranov equation ---
      params = solve_GS(params,1);
      if isempty(params)
        params = [];
        fprintf('       --- FBT did NOT converge\n')
        return
      end
      
      % --- Calculating mutuals between plasma and machine ---
      params = plasma_mutuals(params);
      
      err = nanmin([1e3 sum(abs(Ia./max(Ia)).*abs((params.Ia - Ia)./mean(abs(Ia))))]);
      fprintf('       --- Error: %4.2e\n',err)
      
      if params.iter >= params.Niter, fprintf('   --- Maximum number of params.iterations is %u\n',params.Niter), fprintf('   --- FBT did NOT converge\n'), imsg = 1; end
      
  end
  
  if ~imsg, fprintf('\n  --- FBT converged ---\n'), end

  % --- Writing GEQDSK file ---
   if ~imsg
        params.geqdsk.filename = params.geqdsk.filename(3:end);
          if params.i_save
              if ~isempty(params.shot) && ~isempty(params.eqtime)
                  params.geqdsk.shot     = params.shot;
                  params.geqdsk.time     = params.eqtime/1000;
                  params.geqdsk.filename = sprintf('g%06u.%05u',params.geqdsk.shot,params.geqdsk.time*1000);
                  params.geqdsk.current  = 1/params.geqdsk.rcentr*trapz(linspace(params.geqdsk.simag,params.geqdsk.sibry,params.geqdsk.nh),params.geqdsk.Vp.*params.geqdsk.jtor);
              end
              curdir = pwd;
              fprintf('\n  --- Writing GEQDSK file: %s/%s\n',curdir,params.geqdsk.filename)
              writeg(params.geqdsk,[curdir '/' params.geqdsk.filename])
              fprintf('  --- Writing FBT dataset: %s/p%06u.%05u.mat\n',curdir,params.geqdsk.shot,params.geqdsk.time*1000)
              save(sprintf('%s/p%06u.%05u.mat',curdir,params.geqdsk.shot,params.geqdsk.time*1000),'-struct','params')
              %fprintf('  --- Writing coil currents file: %s/current%06u.%05u.mat\n',curdir,params.geqdsk.shot,params.geqdsk.time*1000)
              Ia = params.Ia.*params.coils(:,5)/1e3;
              Ia = [[zeros(7,1); Ia], zeros(7+length(params.coils(:,1)),1)];
              save(sprintf('%s/current%06u.%05u.dat',curdir,params.geqdsk.shot,params.geqdsk.time*1000),'-ascii','Ia')
          end
   end
  
end
