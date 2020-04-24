function select_point(varargin)
%SELECT_POINT Specify what happens on mouse click
%   check for CBool and add to constraints list
global CBool
global dcm_obj
global C
if CBool == 1
    info_struct = getCursorInfo(dcm_obj);
    C = [C;info_struct.Position];
    disp('Constraints selected:')
    disp(C)
end
end

