function [v] = deformation4(H,Hidx,C,obj)
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
Atop = zeros(2*Fnum, Vnum*2);
G_mat = zeros(4, 6, Fnum);
for i = 1:Fnum
    connected = obj.f(i,:);
    Con1 = [connected(1), connected(2)];
    Con2 = [connected(1), connected(3)];
    Con3 = [connected(2), connected(3)];
    % get neighbours coords first vector of patch
    neigh = [obj.v(connected(1),:);obj.v(connected(2),:);obj.v(connected(3),:)];
    G_k = [neigh(1,1), neigh(1,2),1,0;neigh(1,2),-neigh(1,1),0,1;...
        neigh(2,1), neigh(2,2),1,0;neigh(2,2),-neigh(2,1),0,1;...
        neigh(3,1),neigh(3,2),1,0;neigh(3,2),-neigh(3,1),0,1];
    G_calc = inv(G_k'*G_k)*G_k';
    %save G_calc for later
    G_mat(:,:,i) = G_calc;
    % for edge  of neigh 1 -> 2
    e_kx = neigh(2,1)-neigh(1,1);
    e_ky = neigh(2,2)-neigh(1,2);
    % h matrix
    h = [-1,0,1,0,0,0;0,-1,0,1,0,0]-[e_kx,e_ky;e_ky,-e_kx]...
        *G_calc(1:2,:);
    % Atop
    k = i*2-1;
    Atop(k, Con1(1)*2 -1) = h(1,1);
    Atop(k, Con1(1)*2) = h(1,2);
    Atop(k+1,Con1(1)*2-1) = h(2,1);
    Atop(k+1, Con1(1)*2) = h(2,2);
    Atop(k, Con1(2)*2 -1) = h(1,3);
    Atop(k, Con1(2)*2) = h(1,4);
    Atop(k+1,Con1(2)*2-1) = h(2,3);
    Atop(k+1, Con1(2)*2) = h(2,4);
    Atop(k, Con2(2)*2 -1) = h(1,5);
    Atop(k, Con2(2)*2) = h(1,6);
    Atop(k+1,Con2(2)*2-1) = h(2,5);
    Atop(k+1, Con2(2)*2) = h(2,6);
    
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
eb = zeros(Fnum*2,1);
b = [eb;cb];

%% solve for v
v = (A'*A)\(A'*b);
x=[];
y=[];
for i = 1:size(v,1)
    if mod(i,2) ~= 0 
        x = [x;v(i)];
    else
        y = [y;v(i)];
    end
end
%v = [x,y,zeros(size(x,1),1)];
v = [x,y];
% do scale addjustment
T = zeros(2,2,Fnum);
%use the previously computed G to recalc v' with scale adjust
for i = 1:Fnum
    G = G_mat(:,:,i);
    connected = obj.f(i,:);
    Con1 = [connected(1), connected(2)];
    Con2 = [connected(1), connected(3)];
    Con3 = [connected(2), connected(3)];
    vec = [v(Con1(1),:)';v(Con1(2),:)';v(Con2(2),:)'];
    cs=G(1:2,:)*full(vec);
    c = cs(1);
    s = cs(2);
    T(:,:,i) = (1/(c^2 + s^2)).*[c,s;-s,c];
end
v = scale_correction(H,Hidx,C,obj,T);
end

