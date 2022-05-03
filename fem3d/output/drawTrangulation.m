%clear all
close all

load elements.dat
load nodes.dat

row=size(elements,1);

for ii=1:row
    element=elements(ii,1:4);
    for jj=1:4
        for kk=(jj+1):4
            line(nodes(element([jj,kk]),1),nodes(element([jj,kk]),2),nodes(element([jj,kk]),3),'LineStyle','-')
        end
    end
    vol(ii)=det([[1;1;1;1],nodes(element,1:3)])/6;
end

nr=size(nodes,1);

hold on
index=(nodes(:,4)==-1);
plot3(nodes(index,1),nodes(index,2),nodes(index,3),'ro')


az = 15; el = 30;
view(az,el)
axis square
axis off