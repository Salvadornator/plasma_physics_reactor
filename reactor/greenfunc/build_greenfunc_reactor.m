function build_greenfunc_reactor

fprintf('\n')
fprintf('   --- Start building Green''s functions tables ---')

%% Defining grid size
isave = 1;
Ngrid = [65 129];
for ing = 1:2
    
    Ng = Ngrid(ing);
    % Starting calculation 
    R1 = 0.5;
    R2 = 6;
    dZ = 4;
    Z0 = 0;
    NR = 5;

        %% Clearing figure
        figure(1)
        clf
        plot_reactor
        drawnow

        %% Calculating grid points
        out.r = linspace(R1,R2,Ng);
        out.z = linspace(Z0-dZ,Z0+dZ,Ng);
        [R,Z] = meshgrid(out.r,out.z);

        %% Loading machine coils
        W = load('/home/fmsalvador/matlab/plasma_physics_project/reactor/reactor_vacuum_vessel');
        C = load('/home/fmsalvador/matlab/plasma_physics_project/reactor/reactor_coils');
        cnames = {'OH' 'OH1' 'OH2' 'OH3' 'OH4' 'OH5' 'OH6' 'OH7' 'E1' 'E2' 'E3' 'E4' 'E5' 'E6' 'E7' 'E8' 'E9' 'E10' 'F1' 'F2' 'F3' 'F4' 'D1' 'D2' 'D3' 'D4'};
        out.OH.G = [];

        %% Calculating Green's functions between external coils and grid points
        for icoils = 1:size(C,1)
            Ro = C(icoils,1);
            Zo = C(icoils,2);
            dR = C(icoils,3);
            dZ = C(icoils,4);
            NZ = round(C(icoils,5)/NR);
            Rin = Ro - dR/2 + dR/NR/2;
            Rout = Ro + dR/2 - dR/NR/2;
            Zin = Zo - dZ/2 + dZ/NZ/2;
            Zout = Zo + dZ/2 - dZ/NZ/2;
            Rfilaments = linspace(Rin,Rout,NR);
            Zfilaments = linspace(Zin,Zout,NZ);

            G = zeros(size(R)); dMdR = G; dMdZ = G; BR = G; BZ = G; dBRdR = G; dBZdR = G; dBRdZ = G; dBZdZ = G;
            fprintf('\n')
            fprintf('   --- Coil OH1: [R Z] =  [00 000]')
            for ir = 1:NR
                for iz = 1:NZ
                    fprintf(repmat('\b', 1, 31+length(cnames{icoils+1})))
                    fprintf('   --- Coil %s: [R Z] =  [%2.0u %3.0u]',cnames{icoils+1},ir,iz)
                    G  = G  + C(icoils,5)/NR/NZ*green_em(Rfilaments(ir),Zfilaments(iz),R,Z,'psi');
                    BR = BR + C(icoils,5)/NR/NZ*green_em(Rfilaments(ir),Zfilaments(iz),R,Z,'br');
                    BZ = BZ + C(icoils,5)/NR/NZ*green_em(Rfilaments(ir),Zfilaments(iz),R,Z,'bz');
                    dMdR = dMdR + C(icoils,5)/NR/NZ*green_em(Rfilaments(ir),Zfilaments(iz),R,Z,'dMdr');
                    dMdZ = dMdZ + C(icoils,5)/NR/NZ*green_em(Rfilaments(ir),Zfilaments(iz),R,Z,'dMdz');
                    dBRdR = dBRdR + C(icoils,5)/NR/NZ*green_em(Rfilaments(ir),Zfilaments(iz),R,Z,'dbrdr');
                    dBZdR = dBZdR + C(icoils,5)/NR/NZ*green_em(Rfilaments(ir),Zfilaments(iz),R,Z,'dbzdr');
                    dBRdZ = dBRdZ + C(icoils,5)/NR/NZ*green_em(Rfilaments(ir),Zfilaments(iz),R,Z,'dbrdz');
                    dBZdZ = dBZdZ + C(icoils,5)/NR/NZ*green_em(Rfilaments(ir),Zfilaments(iz),R,Z,'dbzdz');
                end
            end
        eval(['out.' cnames{icoils+1} '.G = G;'])
        eval(['out.' cnames{icoils+1} '.BR = BR;'])
        eval(['out.' cnames{icoils+1} '.BZ = BZ;'])
        eval(['out.' cnames{icoils+1} '.dMdR = dMdR;'])
        eval(['out.' cnames{icoils+1} '.dMdZ = dMdZ;'])
        eval(['out.' cnames{icoils+1} '.dBRdR = dBRdR;'])
        eval(['out.' cnames{icoils+1} '.dBZdR = dBZdR;'])
        eval(['out.' cnames{icoils+1} '.dBRdZ = dBRdZ;'])
        eval(['out.' cnames{icoils+1} '.dBZdZ = dBZdZ;'])

        % Plotting
        figure(1)
        clf
        plot_reactor
        contour(out.r,out.z,G,51,'b')
        for ir = 1:NR
            for iz = 1:NZ
                plot(Rfilaments(ir),Zfilaments(iz),'ok')
            end
        end
        xlabel('R (m)')
        ylabel('Z (m)')
        drawnow
        pause(2)

        end
    out.OH.G = out.OH1.G + out.OH2.G + out.OH3.G + out.OH4.G + out.OH5.G + out.OH6.G + out.OH7.G;
    out.OH.BR = out.OH1.BR + out.OH2.BR + out.OH3.BR + out.OH4.BR + out.OH5.BR + out.OH6.BR + out.OH7.BR;
    out.OH.BZ = out.OH1.BZ + out.OH2.BZ + out.OH3.BZ + out.OH4.BZ + out.OH5.BZ + out.OH6.BZ + out.OH7.BZ;
    out.OH.dMdR = out.OH1.dMdR + out.OH2.dMdR + out.OH3.dMdR + out.OH4.dMdR + out.OH5.dMdR + out.OH6.dMdR + out.OH7.dMdR;
    out.OH.dMdZ = out.OH1.dMdZ + out.OH2.dMdZ + out.OH3.dMdZ + out.OH4.dMdZ + out.OH5.dMdZ + out.OH6.dMdZ + out.OH7.dMdZ;
    out.OH.dBRdR = out.OH1.dBRdR + out.OH2.dBRdR + out.OH3.dBRdR + out.OH4.dBRdR + out.OH5.dBRdR + out.OH6.dBRdR + out.OH7.dBRdR;
    out.OH.dBZdR = out.OH1.dBZdR + out.OH2.dBZdR + out.OH3.dBZdR + out.OH4.dBZdR + out.OH5.dBZdR + out.OH6.dBZdR + out.OH7.dBZdR;
    out.OH.dBRdZ = out.OH1.dBRdZ + out.OH2.dBRdZ + out.OH3.dBRdZ + out.OH4.dBRdZ + out.OH5.dBRdZ + out.OH6.dBRdZ + out.OH7.dBRdZ;
    out.OH.dBZdZ = out.OH1.dBZdZ + out.OH2.dBZdZ + out.OH3.dBZdZ + out.OH4.dBZdZ + out.OH5.dBZdZ + out.OH6.dBZdZ + out.OH7.dBZdZ;

    %% Saving data
    if isave
        fname = ['/home/fmsalvador/matlab/plasma_physics_project/reactor/greenfunc/green_table_' num2str(Ng) 'x' num2str(Ng) '_fbt.mat'];
        fprintf('\n\n')
        fprintf('   --- Saving Green''s functions table saved at %s\n',fname)
        save(fname,'-struct','out')
    end
    fprintf('\n   --- Done ---\n')
end

end