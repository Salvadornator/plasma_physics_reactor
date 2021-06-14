function writeg(gfile,filename)

gfile.psirz  = gfile.psirz.';
gfile.psirzn = gfile.psirzn.';
gfile.sibry  = gfile.sibry;
gfile.simag  = gfile.simag;

fprintfE0 = @(fileid,formatspec,num) fprintf(fileid,'%s',regexprep(sprintf(formatspec,num*10),'(\d)(\.)(.*)','0$2$1$3'));

fid = fopen(filename,'w');
fprintf(fid,'  %-9s',gfile.fittype);
fprintf(fid,'%-14s',gfile.date);
fprintf(fid,'%s','#');
fprintf(fid,'%06d ',gfile.shot);
fprintf(fid,'%05d',gfile.time*1000);
fprintf(fid,'%-13s','ms');
fprintf(fid,'%d',gfile.idum);
fprintf(fid,'%4d',gfile.nw);
fprintf(fid,'%4d\n',gfile.nh);
fprintfE0(fid,'%15.8E',gfile.rdim);
fprintfE0(fid,'%15.8E',gfile.zdim);
fprintfE0(fid,'%15.8E',gfile.rcentr);
fprintfE0(fid,'%15.8E',gfile.rleft);
fprintfE0(fid,'%15.8E\n',gfile.zmid);
fprintfE0(fid,'%15.8E',gfile.rmaxis);
fprintfE0(fid,'%15.8E',gfile.zmaxis);
fprintfE0(fid,'%15.8E',gfile.simag);
fprintfE0(fid,'%15.8E',gfile.sibry);
fprintfE0(fid,'%15.8E\n',gfile.bcentr);
fprintfE0(fid,'%15.8E',gfile.current);
fprintfE0(fid,'%15.8E',gfile.simag);
fprintfE0(fid,'%15.8E',0);
fprintfE0(fid,'%15.8E',gfile.rmaxis);
fprintfE0(fid,'%15.8E\n',0);
fprintfE0(fid,'%15.8E',gfile.zmaxis);
fprintfE0(fid,'%15.8E',0);
fprintfE0(fid,'%15.8E',gfile.sibry);
fprintfE0(fid,'%15.8E',0);
fprintfE0(fid,'%15.8E\n',0);

% Writing fpol
for ii = 1:fix(gfile.nw/5)
    for jj = 1:5
        fprintfE0(fid,'%15.8E',gfile.fpol(5*(ii-1) + jj));
    end
    fprintf(fid,'\n');
end
for kk = 1:rem(gfile.nw,5)
    fprintfE0(fid,'%15.8E',gfile.fpol(fix(gfile.nw/5)*5 + kk));
end
fprintf(fid,'\n');

% Writing pres
for ii = 1:fix(gfile.nw/5)
    for jj = 1:5
        fprintfE0(fid,'%15.8E',gfile.pres(5*(ii-1) + jj));
    end
    fprintf(fid,'\n');
end
for kk = 1:rem(gfile.nw,5)
    fprintfE0(fid,'%15.8E',gfile.pres(fix(gfile.nw/5)*5 + kk));
end
fprintf(fid,'\n');

% Writing ffprim
for ii = 1:fix(gfile.nw/5)
    for jj = 1:5
        fprintfE0(fid,'%15.8E',gfile.ffprim(5*(ii-1) + jj));
    end
    fprintf(fid,'\n');
end
for kk = 1:rem(gfile.nw,5)
    fprintfE0(fid,'%15.8E',gfile.ffprim(fix(gfile.nw/5)*5 + kk));
end
fprintf(fid,'\n');

% Writing pprime
for ii = 1:fix(gfile.nw/5)
    for jj = 1:5
        fprintfE0(fid,'%15.8E',gfile.pprime(5*(ii-1) + jj));
    end
    fprintf(fid,'\n');
end
for kk = 1:rem(gfile.nw,5)
    fprintfE0(fid,'%15.8E',gfile.pprime(fix(gfile.nw/5)*5 + kk));
end
fprintf(fid,'\n');

% Writing psirz
psirz = gfile.psirz(:);
for ii = 1:fix(gfile.nw*gfile.nh/5)
    for jj = 1:5
        fprintfE0(fid,'%15.8E',psirz(5*(ii-1) + jj));
    end
    fprintf(fid,'\n');
end
for kk = 1:rem(gfile.nw*gfile.nh,5)
    fprintfE0(fid,'%15.8E',psirz(fix(gfile.nw*gfile.nh/5)*5 + kk));
end
fprintf(fid,'\n');

% Writing qpsi
for ii = 1:fix(gfile.nw/5)
    for jj = 1:5
        fprintfE0(fid,'%15.8E',gfile.qpsi(5*(ii-1) + jj));
    end
    fprintf(fid,'\n');
end
for kk = 1:rem(gfile.nw,5)
    fprintfE0(fid,'%15.8E',gfile.qpsi(fix(gfile.nw/5)*5 + kk));
end
fprintf(fid,'\n');

% Writting nbbbs and limitr
fprintf(fid,'   %-5d',gfile.nbbbs);
fprintf(fid,'%d',gfile.limitr);
fprintf(fid,'\n');

% Writing rbbbs and zbbbs
ss = 0;
for ii = 1:fix(gfile.nbbbs/5)
    for jj = 1:5
        fprintfE0(fid,'%15.8E',gfile.rbbbs(5*(ii-1) + jj));
        ss = ss + 1;
        if rem(ss,5) == 0
            fprintf(fid,'\n');
        end
        fprintfE0(fid,'%15.8E',gfile.zbbbs(5*(ii-1) + jj));
        ss = ss + 1;
        if rem(ss,5) == 0
            fprintf(fid,'\n');
        end
    end
end
for kk = 1:rem(gfile.nbbbs,5)
    fprintfE0(fid,'%15.8E',gfile.rbbbs(fix(gfile.nbbbs/5)*5 + kk));
    ss = ss + 1;
        if rem(ss,5) == 0
            fprintf(fid,'\n');
        end
    fprintfE0(fid,'%15.8E',gfile.zbbbs(fix(gfile.nbbbs/5)*5 + kk));
    ss = ss + 1;
        if rem(ss,5) == 0
            fprintf(fid,'\n');
        end
end
fprintf(fid,'\n');

% Writing rlim and zlim
ss = 0;
for ii = 1:fix(gfile.limitr/5)
    for jj = 1:5
        fprintfE0(fid,'%15.8E',gfile.rlim(5*(ii-1) + jj));
        ss = ss + 1;
        if rem(ss,5) == 0
            fprintf(fid,'\n');
        end
        fprintfE0(fid,'%15.8E',gfile.zlim(5*(ii-1) + jj));
        ss = ss + 1;
        if rem(ss,5) == 0
            fprintf(fid,'\n');
        end
    end
end
for kk = 1:rem(gfile.limitr,5)
    fprintfE0(fid,'%15.8E',gfile.rlim(fix(gfile.limitr/5)*5 + kk));
    ss = ss + 1;
        if rem(ss,5) == 0
            fprintf(fid,'\n');
        end
    fprintfE0(fid,'%15.8E',gfile.zlim(fix(gfile.limitr/5)*5 + kk));
    ss = ss + 1;
        if rem(ss,5) == 0
            fprintf(fid,'\n');
        end
end
fprintf(fid,'\n');

% Writing kvtor, rvtor and nmass
fprintf(fid,'  %3d',gfile.kvtor);
fprintfE0(fid,'%15.8E',gfile.rvtor);
fprintf(fid,'%5d\n',gfile.nmass);
    
if gfile.kvtor > 0

    % Writing pressw
    for ii = 1:fix(gfile.nw/5)
        for jj = 1:5
            fprintfE0(fid,'%15.8E',gfile.pressw(5*(ii-1) + jj));
        end
        fprintf(fid,'\n');
    end
    for kk = 1:rem(gfile.nw,5)
        fprintfE0(fid,'%15.8E',gfile.pressw(fix(gfile.nw/5)*5 + kk));
    end
    fprintf(fid,'\n');

    % Writing pwprim
    for ii = 1:fix(gfile.nw/5)
        for jj = 1:5
            fprintfE0(fid,'%15.8E',gfile.pwprim(5*(ii-1) + jj));
        end
        fprintf(fid,'\n');
    end
    for kk = 1:rem(gfile.nw,5)
        fprintfE0(fid,'%15.8E',gfile.pwprim(fix(gfile.nw/5)*5 + kk));
    end
    fprintf(fid,'\n');
end

if gfile.nmass > 0
    % Writing dmion
    for ii = 1:fix(gfile.nw/5)
        for jj = 1:5
            fprintfE0(fid,'%15.8E',gfile.dmion(5*(ii-1) + jj));
        end
        fprintf(fid,'\n');
    end
    for kk = 1:rem(gfile.nw,5)
        fprintfE0(fid,'%15.8E',gfile.dmion(fix(gfile.nw/5)*5 + kk));
    end
    fprintf(fid,'\n');
end

% Writing rhovn
if ~isfield(gfile,'rhovn')
    rhov = integration(linspace(0,1,gfile.nw),gfile.qpsi',0);
    gfile.rhovn = sqrt((rhov - rhov(1))/(rhov(end) - rhov(1)));
end
for ii = 1:fix(gfile.nw/5)
    for jj = 1:5
        fprintfE0(fid,'%15.8E',gfile.rhovn(5*(ii-1) + jj));
    end
    fprintf(fid,'\n');
end
for kk = 1:rem(gfile.nw,5)
    fprintfE0(fid,'%15.8E',gfile.rhovn(fix(gfile.nw/5)*5 + kk));
end
fprintf(fid,'\n');
fprintf(fid,'    %d',0);
fprintf(fid,'\n');
fclose(fid);
%disp(['Writing GEQDSK file: ' filename])

end