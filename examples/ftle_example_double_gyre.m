% Example use of Coherent_structure_starter_pack by showing computation of
% FTLE for the double gyre model
clear; clc;
addpath('../functions')
%% Define xygrid, parameters and velocity field
set(0,'defaultAxesFontSize',24);
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaulttextInterpreter','latex');
make_movie = false; % change to true to produce and save a movie of FTLE
nx = 201;
ny = 101;
x = linspace(0,2,nx);
y = linspace(0,1,ny);
[X,Y] = meshgrid(x,y);
t0 = 0;
T = 10;
intspan = [t0,T];
% Make tspan longer to produce a series of FTLE data/plots, currently only
% setup to display the FTLE field at time t = t0, required to make movie
tspan = linspace(t0,T,501);
vfield = @double_gyre;
ftle = FTLE_field_fast(vfield,tspan,x,y,T);
%% Produce plots of the ftle field for the double gyre
figure
for k = 1:size(ftle,3)
    ncontourf(X,Y,ftle(:,:,1),40)
    title(['Double Grye FTLE at t = ',num2str(tspan(k))]);
    xlabel('x');
    ylabel('y');
    mov(k) = getframe(gcf);
end
%% Make movie if make_movie set to true
if make_movie
    fps = 10;
    obj = VideoWriter('double_gyre_ftle_T10.avi');
    obj.FrameRate = fps;
    open(obj);
    writeVideo(obj,mov);
    close(obj);
end 