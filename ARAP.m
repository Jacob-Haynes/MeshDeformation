%% This script runs the ARAP algorithm
% put info here

% Clear
clc; clear all;
% set global vars
global CBool
global C
global dcm_obj
global HBool
global H
global Hclick
global obj
global Vnum
global Fnum
global arap
%% Read mesh
obj = readmesh('man.obj');
% number of verticies and faces
Vnum = size(obj.v,1);
Fnum = size(obj.f,1);

%% display model
% create global bool for selecting constraints
CBool = 0;
HBool = 0;
Hclick = 0;
% set constraints global
C = [];
% create gui
%create figure
fig = figure('visible','off');
ax = axes('units','pixels');
hold on
set(ax, 'Yticklabel',[]);
set(ax, 'Xticklabel',[]);
%create togle select constraints
global toggleConst
toggleConst = uicontrol('Style', 'togglebutton', 'String',...
    'Chose Constraints', 'position' , [20 20 100 20],...
    'Callback', @tc);
fig.Visible = 'on';
%display the mesh
dispmodel(obj);

% select point button
global togglePos
togglePos = uicontrol('Style', 'pushbutton', 'String',...
    'Select Point', 'position' , [150 20 100 20],...
    'Callback', @select_point);

% Select deformation handle
global toggleHand
toggleHand = uicontrol('Style', 'togglebutton', 'String',...
    'Move Handle', 'position' , [280 20 100 20],...
    'Callback', @move_handle);

% click and move functionality
handles = guidata(gca);
set(gcf, 'windowbuttondownfcn', {@myclick});

