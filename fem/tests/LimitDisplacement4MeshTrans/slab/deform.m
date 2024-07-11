load mesh.nodes;
f=mesh;
R=2.0;
Z=1.0;
f(:,5)=Z-sqrt(max(0,R^2-f(:,3).^2));
fid=fopen("mesh.nodes2","w")
for i=1:size(f,1)
	fprintf(fid,"%d -1 %12.6e %12.6e %12.6e\n",f(i,1),f(i,3),f(i,4),f(i,5))
endfor
fclose(fid);
