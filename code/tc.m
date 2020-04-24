function tc(toggleConst, eventData)
%TC when constraints button is pressed get user to supply a constraint
%point
global CBool
global dcm_obj
global C

if CBool == 0
    CBool = 1;
    disp('Click constraint points, click this again to stop selecting constraints')
    disp('Click SELECT POINT button to confirm a point as constraint')
    datacursormode on
    dcm_obj = datacursormode(gcf);
    set(dcm_obj, 'DisplayStyle', 'Window') %change to window
else
    CBool = 0;
    disp('Done selecting constraint')
    datacursormode off
    disp('Constraints Selected:')
    disp(C)
    %plot the constraints
    scatter(gca, C(:,1), C(:,2),'r', 'filled')
end

end

