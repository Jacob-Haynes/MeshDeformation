function dispmodel(obj)
%DISPMODEL displays obj struct
% display object
patch('vertices', obj.v, 'faces', obj.f,'EdgeColor',[1 0.4 0.6],'FaceColor','none','LineWidth',1);
end

