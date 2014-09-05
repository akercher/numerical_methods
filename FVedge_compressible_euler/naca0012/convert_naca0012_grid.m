% convert_naca0012_grid.m
clear all;

fid = fopen('naca0012.mesh.coarse.backup');
% $$$ fid = fopen('naca0012.mesh.medium.backup');
% $$$ fid = fopen('naca0012.mesh.fine.backup');

nnode = fscanf(fid,'%f',1);

cbuffer = fscanf(fid,'%s',1);
cbuffer = fscanf(fid,'%s',[1,3]);
cbuffer = fscanf(fid,'%s',1);
cbuffer = fscanf(fid,'%s',[1,2]);

% $$$ array = fscanf(fid,'%f',1)
ndim = fscanf(fid,'%f',1);
ntype = fscanf(fid,'%f',1);

cbuffer = fscanf(fid,'%s',[1,4]);
nelem = fscanf(fid,'%f',1);
npoin = fscanf(fid,'%f',1);
nboun = fscanf(fid,'%f',1);
fbuffer = fscanf(fid,'%f',1);

%read in connectivity
cbuffer = fscanf(fid,'%s',[1,7]);
for ielem=1:nelem
  ibuffer = fscanf(fid,'%f',1);
  lnode(ielem,:) = fscanf(fid,'%f',[1,nnode]);
  fbuffer = fscanf(fid,'%f',[1,8]);
end

%read in point coordinates
cbuffer = fscanf(fid,'%s',[1,4]);
for ipoin=1:npoin
  ibuffer = fscanf(fid,'%f',1);
  lpoin(ipoin,:) = fscanf(fid,'%f',[1,ndim]);
end

%read initial valuees
cbuffer = fscanf(fid,'%s',[1,8]);
for ipoin=1:npoin
  fbuffer = fscanf(fid,'%f',[1,12]);  
end

%read in boundary conditions
cbuffer = fscanf(fid,'%s',[1,2]);

%external and airfoil boundary nodes
bext = [];
bwall = [];
for iboun=1:nboun
  bcond(iboun,:) = fscanf(fid,'%f',[1,2]);
  fbuffer = fscanf(fid,'%f',[1,5]);
  if bcond(iboun,2) == 4
    bext = [bext;bcond(iboun,1)];
  else  
    bwall = [bwall;bcond(iboun,1)];    
  end    
end  
fclose(fid);

%boundaries are closed --> add first bnode to end of list
bext = [bext;bext(1)];
bwall = [bwall;bwall(1)];

fid = fopen('naca0012.mesh.coarse','w');
% $$$ fid = fopen('naca0012.mesh.medium','w');
% $$$ fid = fopen('naca0012.mesh.fine','w');

% number of boundary conditions
nbcs = 2;
nbnode = [length(bext),length(bwall)];

fprintf(fid,'%d %d %d\n',npoin,nelem,0);

%write coordinates
for i=1:npoin
  fprintf(fid,'%f %f\n',lpoin(i,1),lpoin(i,2));
end

%write connectivity
for i=1:nelem
  fprintf(fid,'%d %d %d\n',lnode(i,1),lnode(i,2),lnode(i,3));
end

%write number of boundary conditions and 
%number of boundary nodes for each condition
fprintf(fid,'%d\n',nbcs);
for i=1:nbcs
  fprintf(fid,'%d\n',nbnode(i));
end

%write boundary nodes
for i=1:nbnode(1)
  fprintf(fid,'%d\n',bext(i));
end
for i=1:nbnode(2)
  fprintf(fid,'%d\n',bwall(i));
end

