function myclick(h,event,type)
global HBool
global Hclick
global obj
global H
global C
global Hidx
global arap
if HBool == 1
    if Hclick == 0
        %select point to move thats nearest cursor click
        out = get(gca,'CurrentPoint'); %get mouse click
        Xcoord = obj.v(:,1); % get coords of mesh
        Ycoord = obj.v(:,2);
        aspect = daspect(gca); %aspect ratio of axis
        dist = sqrt(((out(1,1)-Xcoord).*aspect(2)).^2+((out(1,2)-Ycoord).*aspect(1)).^2); %dist of click to each mesh point
        idx=find(dist == min(dist),1,'first'); %get index of closest point
        H = obj.v(idx,:);
        %then set Hclick to 1
        Hclick = 1;
        % give user instructions
        disp('Point selected:')
        disp(H)
        disp('Now choose the destination')
    else
        %chose new mouse click position 
        out = get(gca,'CurrentPoint');
        out = out(1,:); %z to match the mesh
        out(1,3)=0; 
        disp('Target Destination:')
        disp(out)
        %set the handle new position to click position
        Hidx = find(obj.v ==H,1,'first');
        H=out;
        %then set Hclick to 0
        Hclick = 0;
        % perform deformation
        arap.v = deformation4(H,Hidx,C,obj);
        arap.v = full(arap.v);
        arap.f=obj.f;
        cla
        dispmodel(arap);
    end
end
end


