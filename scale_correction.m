function [v] = scale_correction(H,Hidx,C,obj,T)
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
Atop = zeros(Vnum^2, Vnum);
e = zeros(2,1,Vnum^2);
e_x = zeros(Vnum^2,1); %also do edge matrix in same loop for speed
e_y = zeros(size(e_x));
for i = 1:Fnum
    connected = obj.f(i,:);
    Con1 = [connected(1), connected(2)];
    Con2 = [connected(1), connected(3)];
    Con3 = [connected(2), connected(3)];
    Atop(Con1(1)*Con1(2), Con1(1)) = -1;
    Atop(Con1(1)*Con1(2), Con1(2)) = 1;
    Atop(Con2(1)*Con2(2), Con2(1)) = -1;
    Atop(Con2(1)*Con2(2), Con2(2)) = 1;
    Atop(Con3(1)*Con3(2), Con3(1)) = -1;
    Atop(Con3(1)*Con3(2), Con3(2)) = 1;
    %do edge matrix in same loop for speed 
    %include T
    e_x(Con1(1)*Con1(2)) = (obj.v(Con1(2),1)-obj.v(Con1(1),1));
    e_x(Con2(1)*Con2(2)) = (obj.v(Con2(2),1)-obj.v(Con2(1),1));
    e_x(Con3(1)*Con3(2)) = (obj.v(Con3(2),1)-obj.v(Con3(1),1));
    e_y(Con1(1)*Con1(2)) = (obj.v(Con1(2),2)-obj.v(Con1(1),2));
    e_y(Con2(1)*Con2(2)) = (obj.v(Con2(2),2)-obj.v(Con2(1),2));
    e_y(Con3(1)*Con3(2)) = (obj.v(Con3(2),2)-obj.v(Con3(1),2));
    
    e(:,:,Con1(1)*Con1(2)) = T(:,:,i)*[e_x(Con1(1)*Con1(2)), e_y(Con1(1)*Con1(2))]';
    e(:,:,Con2(1)*Con2(2)) = T(:,:,i)*[e_x(Con2(1)*Con2(2)), e_y(Con2(1)*Con2(2))]';
    e(:,:,Con3(1)*Con3(2)) = T(:,:,i)*[e_x(Con3(1)*Con3(2)), e_y(Con3(1)*Con3(2))]';
end
for i = 1:Vnum^2
    e_x(i) = e(1,1,i);
    e_y(i) = e(2,1,i);
end
Atop = sparse(Atop);
% build bottom of A
Abot = zeros(Vnum);
for i = 1:size(Cidx,1)
    Abot(Cidx(i),Cidx(i)) = w;
end
Abot = sparse(Abot);
%build A
A = [Atop;Abot];

%% construct b
% build constraints matrix
c_x = zeros(Vnum,1);
c_y = c_x;
c_x(Cidx) = C(:,1);
c_y(Cidx) = C(:,2);
c_x(Hidx) = H(:,1);
c_y(Hidx) = H(:,2);
c_x = w.*c_x;
c_y = w.*c_y;
c_x = sparse(c_x);
c_y = sparse(c_y);

%combine edge and constraints
b_x = [e_x;c_x];
b_y = [e_y;c_y];


%% solve for v
v_x = (A'*A)\(A'*b_x);
v_y = (A'*A)\(A'*b_y);
v = [v_x,v_y,zeros(size(v_x))];

end

