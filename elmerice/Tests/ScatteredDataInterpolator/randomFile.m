%% matlab function to generate sparse datasets
%%% 1 with random point positions
%%% 1 along 2 fictive "flight" lines
function randomFile

%%%%% generate a data set with random point positions
xx=rand(200,2);

xyz(:,1)=210e03*xx(:,1)-5e03;
xyz(:,2)=60e03*xx(:,2)-5e03;

for ii=1:size(xx,1)
 xyz(ii,3)=zs(xyz(ii,1),xyz(ii,2));
end

save('Rand200.txt','xyz','-ASCII');

%%%%% generate a data set along two "flight" lines at y=15km and y=30km
xgrid=-1000:500:201000;
n1=size(xgrid,2);

xyz2(1:n1,1)=xgrid;
xyz2(1:n1,2)=15000;


xyz2(n1+1:n1+n1,1)=xgrid;
xyz2(n1+1:n1+n1,2)=30000;

for ii=1:size(xyz2,1)
 xyz2(ii,3)=zs(xyz2(ii,1),xyz2(ii,2));
end

save('FlightLines.txt','xyz2','-ASCII');


%%%% The "True" variable
function zs = zs(x,y) 
  Lx = 200.0e3;
  Ly = 50.0e03;
  zs=500.0-1.0e-03*x+20.0*(sin(3.0*pi*x/Lx)*sin(2.0*pi*y/Ly));
end

end
