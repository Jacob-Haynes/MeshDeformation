function obj = readmesh(filename)
%READMESH read obj file to return matlab obj
%   reads verticies and faces
%   INPUT: file name
%   Output: v verticies
%           f faces (3 verticies assumed)

% fields
v = []; f = [];
% load
fileID = fopen(filename);
% go through obj file
while 1
    line = fgetl(fileID);
    %check if end of file
    if ~ischar(line), break, end
    %type of line
    lt = sscanf(line,'%s',1); %read the line type from fist char of line
    switch lt
        case 'v' %if mesh vertex
            v = [v;sscanf(line(2:end),'%f')']; %store the vertex coord 
        case 'f' %if a face
            f = [f;sscanf(line(2:end),'%f')']; %store face
    end
end
fclose(fileID);
% build matlab obj
obj.v = v; obj.f = f;
end