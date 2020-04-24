function move_handle(varargin)
%MOVE_HANDLE Summary of this function goes here
%   Detailed explanation goes here
global HBool
if HBool == 0
    HBool = 1;
    disp('Click and drag deformation handle')

else
    HBool = 0;
    disp('Done deformation')

end
end

