clc; clear; close all;

inp_file = "viewfactor_result.dat";
out_file = "Radiation_InitVF_Prefs.txt";

%% ------------------------------------------------------------------------
% 公式1
d = 2:0.1:8;
vf = zeros(1,length(d));

fun1 = @(c, d) (2.*c - sin(2.*c)) ./ sqrt(d.*d - 4.*cos(c).*cos(c)) .* sin(2.*c);
fun2 = @(c, d) sqrt(1 - 4/d/d .* cos(c).^2) .* sin(c).^2;

for i = 1:length(d)
    if d(i) < 2
        q = integral(@(ita) fun1(ita, d(i)), acos(d(i)^2/4), pi/2);
        vf(i) = 4/(2+d(i)) * q / (d(i)*pi) + d(i)^2/16*(d(i)-2);
    else
        q = integral(@(ita) fun1(ita, d(i)), 0, pi/2);
        vf(i) = q / (d(i)*pi);
    end
end


figure
plot(d, vf, LineWidth=2)
xlabel 'd'
ylabel 'vf'
xlim([1.9 5])
grid on
hold on


%% ------------------------------------------------------------------------
% 公式2
d = 2:0.1:4;
vf = zeros(1,length(d));

for i = 1:length(d)
    q = integral(@(c) fun2(c, d(i)), 0, pi/2);
    vf(i) = 0.5 - 2/pi*q;
end

plot(d, vf, LineWidth=2)
hold on


%% ------------------------------------------------------------------------
data = importdata(inp_file);
data = data.data(data.data(:,2) > 0, :);

% vf_voro = data(:, 5);
vf_voro = data(:, 4);
d_voro  = data(:,3)./0.03;

vf_max = 0.5 - 2/pi*integral(@(c) fun2(c, 2), 0, pi/2);

for i=1:length(d_voro)
    if (d_voro(i) <= 2.0)
        vf_voro(i) = min(vf_voro(i), vf_max);
    elseif (d_voro(i) >= 3)
        vf_voro(i) = 0;
    else
        q = integral(@(c) fun2(c, d_voro(i)), 0, pi/2);
        tmp_vf = 0.5 - 2/pi*q;

        if (vf_voro(i) > tmp_vf)
            vf_voro(i) = 0;
        end
%         vf_voro(i) = min(vf_voro(i), tmp_vf);
    end
end


scatter(d_voro, vf_voro, Marker="+", MarkerEdgeColor="Black")


fileID = fopen(out_file, "w");
fprintf(fileID, "%d\n", length(d_voro));
for i=1:length(d_voro)
%     disp(d_voro);
    if (data(i,1) < data(i,2))
        fprintf(fileID, "%d %d %.18f %.18f\n", data(i, 1), data(i, 2), vf_voro(i), d_voro(i));
%         fprintf(fileID, "%d %d %.18f\n", data(i, 1), data(i, 2), vf_voro(i));
    else
        fprintf(fileID, "%d %d %.18f %.18f\n", data(i, 2), data(i, 1), vf_voro(i), d_voro(i));
%         fprintf(fileID, "%d %d %.18f\n", data(i, 2), data(i, 1), vf_voro(i));
    end
end
fclose(fileID);


%--------------------------
figure
histogram(vf_voro)
title('void fraction')

%--------------------------
disp(['max distance: ' num2str(max(d_voro))]);

figure
histogram(d_voro)
title('distance')


