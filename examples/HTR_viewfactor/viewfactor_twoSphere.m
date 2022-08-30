clc; clear; close all;


%% ------------------------------------------------------------------------
d = 1:0.1:8;
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
xlim([1.5 5])
grid on
hold on


%% ------------------------------------------------------------------------
d = 2:0.1:4;
vf = zeros(1,length(d));

for i = 1:length(d)
    q = integral(@(c) fun2(c, d(i)), 0, pi/2);
    vf(i) = 0.5 - 2/pi*q;
end

plot(d, vf, LineWidth=2)
hold on


%% ------------------------------------------------------------------------
data = importdata("viewfactor_result.dat");
data = data.data(data.data(:,2) > 0, :);

% vf_voro = data(:, 5);
vf_voro = data(:, 4);
d_voro  = data(:,3)./0.03;

vf_max = 0.5 - 2/pi*integral(@(c) fun2(c, 2), 0, pi/2);

for i=1:length(d_voro)
    if (d_voro(i) <= 2.0)
        vf_voro(i) = min(vf_voro(i), vf_max);
    else
        q = integral(@(c) fun2(c, d_voro(i)), 0, pi/2);
        tmp_vf = 0.5 - 2/pi*q;

        vf_voro(i) = min(vf_voro(i), tmp_vf);
    end
end


scatter(d_voro, vf_voro, Marker="+", MarkerEdgeColor="Black")


fileID = fopen("viewfactor_result_optimized.dat", "w");
for i=1:length(d_voro)
    fprintf(fileID, "%d %d %.18f %.18f\n", data(i, 1), data(i, 2), d_voro(i), vf_voro(i));
end
fclose(fileID);

figure
histogram(vf_voro)
