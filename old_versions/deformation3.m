function [v] = deformation3(H,Hidx,C,obj)
%DEFORMATION perform the ARAP deformation edge version
%   perform deformation of mesh obj, with constraints C
%   and handle H - Hidx is the index of H in obj.
%   Return the deformed mesh, arap

global Vnum
global Fnum

% get constraints index in obj.v
Cidx = [];
for i = 1:size(C,1)
    Cidx = [Cidx, find(obj.v ==C(i,:),1,'first')];
end
Cidx = Cidx';
w = 1000; %weight factor

%% construct A
% build top of A
eb = zeros(Fnum*3*2,1);
Atop = zeros(Fnum*3*2 , Vnum*2);
for i = 1:Fnum
    connected = obj.f(i,:);
    Edge1 = [connected(1), connected(2)]; 
    Edge2 = [connected(1), connected(3)];
    Edge3 = [connected(2), connected(3)];
    % get neighbours coords of patch
    neigh = [obj.v(connected(1),:);obj.v(connected(2),:);obj.v(connected(3),:)];
    %make G for each edge
    G_1 = [neigh(1,1), neigh(1,2),1,0;neigh(1,2),-neigh(1,1),0,1;...
        neigh(2,1), neigh(2,2),1,0;neigh(2,2),-neigh(2,1),0,1;...
        neigh(3,1),neigh(3,2),1,0;neigh(3,2),-neigh(3,1),0,1]; 
    G_1Calc = inv(G_1'*G_1)*G_1';
    G_2 = [neigh(1,1), neigh(1,2),1,0;neigh(1,2),-neigh(1,1),0,1;...
        neigh(3,1),neigh(3,2),1,0;neigh(3,2),-neigh(3,1),0,1;...
        neigh(2,1), neigh(2,2),1,0;neigh(2,2),-neigh(2,1),0,1];
    G_2Calc = inv(G_2'*G_2)*G_2';
    G_3 = [neigh(2,1), neigh(2,2),1,0;neigh(2,2),-neigh(2,1),0,1;...
        neigh(3,1),neigh(3,2),1,0;neigh(3,2),-neigh(3,1),0,1;...
        neigh(1,1), neigh(1,2),1,0;neigh(1,2),-neigh(1,1),0,1];
    G_3Calc = inv(G_3'*G_3)*G_3';
    % for edge vectors
    e_kx1 = neigh(2,1)-neigh(1,1);
    e_ky1 = neigh(2,2)-neigh(1,2);
    e_kx2 = neigh(3,1)-neigh(1,1);
    e_ky2 = neigh(3,2)-neigh(1,2);
    e_kx3 = neigh(3,1)-neigh(2,1);
    e_ky3 = neigh(3,2)-neigh(2,2);
    % h matrix
    h1 = [-1,0,1,0,0,0;0,-1,0,1,0,0]-[e_kx1,e_ky1;e_ky1,-e_kx1]...
        *G_1Calc(1:2,:);
    h2 = [-1,0,1,0,0,0;0,-1,0,1,0,0]-[e_kx2,e_ky2;e_ky2,-e_kx2]...
        *G_2Calc(1:2,:); 
    h3 = [-1,0,1,0,0,0;0,-1,0,1,0,0]-[e_kx3,e_ky3;e_ky3,-e_kx3]...
        *G_3Calc(1:2,:);
    % Atop
    for j = 1:3 %for each face do each edge
        k = i*3-3+j;
        a = Edge1(1);
        b = Edge1(2);
        c = Edge2(2);
        if j == 1
            h = h1;
        elseif j == 2
            h = h2;
        else
            h = h3;
        end
        Atop(k, a*2 -1) = h(1,1);
        Atop(k, a*2) = h(1,2);
        Atop(k+1,a*2-1) = h(2,1);
        Atop(k+1, a*2) = h(2,2);
        Atop(k, b*2 -1) = h(1,3);
        Atop(k, b*2) = h(1,4);
        Atop(k+1,b*2-1) = h(2,3);
        Atop(k+1, b*2) = h(2,4);
        Atop(k, c*2 -1) = h(1,5);
        Atop(k, c*2) = h(1,6);
        Atop(k+1,c*2-1) = h(2,5);
        Atop(k+1, c*2) = h(2,6);
    end
end
Atop = sparse(Atop);
% build bottom of A
Abot = zeros(Vnum*2);
for i = 1:size(Cidx,1)
    Abot(Cidx(i)*2-1,Cidx(i)*2-1) = w;
    Abot(Cidx(i)*2, Cidx(i)*2) = w;
end
Abot(Hidx*2-1,Hidx*2-1)=w;
Abot(Hidx*2,Hidx*2)=w;
Abot = sparse(Abot);
%build A
A = [Atop;Abot];

%% construct b
% build constraints matrix
cb = zeros(size(obj.v,1)*2,1);
cb(Cidx*2 -1) = C(:,1);
cb(Cidx*2) = C(:,2);
cb(Hidx*2-1) = H(:,1);
cb(Hidx*2) = H(:,2);
cb = w.*cb;
cb = sparse(cb);
%combine edge and constraints
b = [eb;cb];

%% solve for v
v = (A'*A)\(A'*b);
x=[];
y=[];
for i = 1:size(v,1) %reshape v to obj structure
    if mod(i,2) ~= 0 
        x = [x;v(i)];
    else
        y = [y;v(i)];
    end
end
v = [x,y,zeros(size(x,1),1)];
end

