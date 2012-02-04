function write_targets(filename,stats)

xy = [stats.Centroid]';
x = xy(1:2:end);
y = xy(2:2:end);
n   = [stats.Area]';
sumg = round([stats.sumg]');
nx  = round([stats.MajorAxisLength]');
ny  = round([stats.MinorAxisLength]');
np = length(n);
pnr = 0:np-1;
tnr = -1*ones(np,1);

fid = fopen(filename,'w');
fprintf(fid,'%d\n',np);
fprintf(fid,'%4d %9.4f %9.4f %5d %5d %5d %5d %5d\n',[pnr',x,y,n,nx,ny,sumg,tnr]');
fclose(fid);
