clear all
fid = fopen('mesh.dat', 'rb');

NumNod= fscanf(fid, '%d', 1);
Ncoor = fscanf(fid, '%d', 1);
nodes=[];
for i = 1:NumNod
    for j = 1:Ncoor+1
        nodes(i,j)=fscanf(fid, '%f', 1);
    end
end



NT= fscanf(fid, '%d', 1);
Ncoor = fscanf(fid, '%d', 1);
elements=[];
for i = 1:NT
    for j = 1:Ncoor+1
        elements(i,j)=fscanf(fid, '%d', 1);
    end
end

fclose(fid)

for ii=1:NT
    element=elements(ii,1:4);
    for jj=1:4
        for kk=(jj+1):4
            line(nodes(element([jj,kk]),1),nodes(element([jj,kk]),2),nodes(element([jj,kk]),3),'LineStyle','-','Color','black')
        end
    end
    vol(ii)=det([[1;1;1;1],nodes(element,1:3)])/6;
end

nr=size(nodes,1);

%hold on
%index=(nodes(:,4)==-1);
%plot3(nodes(index,1),nodes(index,2),nodes(index,3),'ro')


az = 15; el = 30;
view(az,el)
axis square
axis off