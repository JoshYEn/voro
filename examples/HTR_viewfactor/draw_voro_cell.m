clc; clear; close all;

%% -------------------------------------------------
vertices = importdata("temp_voronoi_vertices.dat");
% facevtid = importdata("temp_voronoi_facevtid.dat");
fileID = fopen("temp_voronoi_facevtid.dat", "r");
faceIDs = fscanf(fileID, '%d');


%% -------------------------------------------------
figure()
[s_x, s_y, s_z] = sphere(50);
s_x = s_x * 0.03 + vertices(end,2);
s_y = s_y * 0.03 + vertices(end,3);
s_z = s_z * 0.03 + vertices(end,4);

% scatter3(vertices(end,2), vertices(end,3), vertices(end,4), 'filled')
surf = surf(s_x, s_y, s_z, zeros(51));
shading interp
alpha(0.3)
axis equal
view(45,45)
grid on
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Voronoi Cell')


%% -------------------------------------------------
vt_count  = 1;
vt_number = 0;
vt_arr = [];

vt_x = [];
vt_y = [];
vt_z = [];

for i=1:length(faceIDs)
    if (vt_count == 1)
        vt_number = faceIDs(i);
        vt_count = vt_count + 1;
        continue
    end
    if (vt_count < vt_number+1)
        vt_arr = [vt_arr, faceIDs(i)];
        vt_x = [vt_x, vertices(faceIDs(i)+1, 2)];
        vt_y = [vt_y, vertices(faceIDs(i)+1, 3)];
        vt_z = [vt_z, vertices(faceIDs(i)+1, 4)];
        vt_count = vt_count + 1;
    else
        vt_arr = [vt_arr, faceIDs(i)];
        disp(['drawing face:' num2str(vt_arr)]);
        vt_x = [vt_x, vertices(faceIDs(i)+1, 2)];
        vt_y = [vt_y, vertices(faceIDs(i)+1, 3)];
        vt_z = [vt_z, vertices(faceIDs(i)+1, 4)];

        p = patch(vt_x, vt_y, vt_z, "black");
        p.LineWidth = 2;
        alpha(p, 0);
        hold on;
        drawnow;
%         pause;

        vt_arr = [];
        vt_x = []; vt_y = []; vt_z = [];
        vt_count = 1;
    end
end

