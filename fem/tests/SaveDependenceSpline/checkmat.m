load dep.dat
load mat.dat

semilogy(dep(:,2),dep(:,3:5),mat(:,1),mat(:,2),'o')
