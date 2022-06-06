clear all
fid = fopen('mesh.dat', 'rb');

NumNod= fscanf(fid, '%d', 1);
Ncoor = fscanf(fid, '%d', 1);
p=[];
for i = 1:NumNod
    for j = 1:Ncoor
        p(i,j)=fscanf(fid, '%f', 1);
    end
end



NumNod= fscanf(fid, '%d', 1);
Ncoor = fscanf(fid, '%d', 1);
T=[];
for i = 1:NumNod
    for j = 1:Ncoor
        T(i,j)=fscanf(fid, '%d', 1);
    end
end

fclose(fid)

showmesh(p,T)