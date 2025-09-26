clc;clear;

%%
M = 128;
N = 128;
data_name = "imR_RGB_3096_half2full_128x128_6_f = 0.5_1224-082913";
f = 0.5;

%% order
m = 128;
n = 128;
num = floor(m*n*0.1);
y_center = floor((m-1)/2)+1;
x_center = floor((n-1)/2)+1;

[X, Y] = meshgrid(0:n-1, 0:m-1);
values = (X-x_center).^2+(Y-y_center).^2;
values(66:end,:) = values(66:end,:)+m*n;
values(65,1:64) = values(65,1:64)+m*n;
values = ifftshift(values);
[~, I] = sort(values(:));
[row_s, col_s] = ind2sub(size(values), I);

fy_1 = (row_s-1)/m;
fx_1 = (col_s-1)/n;


[X, Y] = meshgrid(0:n-1, 0:m-1);
values = (X-x_center).^2+(Y-y_center).^2;
values(1:64,:) = values(1:64,:)+m*n;
values(65,66:end) = values(65,66:end)+m*n;
values = ifftshift(values);
[~, I] = sort(values(:));
[row_s, col_s] = ind2sub(size(values), I);

fy_2 = (row_s-1)/m;
fx_2 = (col_s-1)/n;

fy = [fy_1(1:num);fy_2(1:num)];
fx = [fx_1(1:num);fx_2(1:num)];


%% Data Processing
data = importdata(strcat(data_name,".txt"));

y = (2*data(1:9:end,:)-data(2:9:end,:)-data(3:9:end,:))...
        +1i*3^0.5*(data(2:9:end,:)-data(3:9:end,:));
x = (2*data(4:9:end,:)-data(5:9:end,:)-data(6:9:end,:))...
        +1i*3^0.5*(data(5:9:end,:)-data(6:9:end,:));
im = (2*data(7:9:end,:)-data(8:9:end,:)-data(9:9:end,:))...
        +1i*3^0.5*(data(8:9:end,:)-data(9:9:end,:));
k = data(1:9:end,:)+data(2:9:end,:)+data(3:9:end,:)...
    +data(4:9:end,:)+data(5:9:end,:)+data(6:9:end,:)...
    +data(7:9:end,:)+data(8:9:end,:)+data(9:9:end,:);


%% locate
y = -angle(y);
y(y<0) = y(y<0)+2*pi;
y = y/2/pi*M/f;

x = -angle(x);
x(x<0) = x(x<0)+2*pi;
x = x/2/pi*N/f;
%% 高斯平滑
window = 5;  % 设置高斯平滑的标准差
x(:,1) = smoothdata(x(:,1),"movmean",window);
x(:,2) = smoothdata(x(:,2),"movmean",window);
x(:,3) = smoothdata(x(:,3),"movmean",window);
y(:,1) = smoothdata(y(:,1),"movmean",window);
y(:,2) = smoothdata(y(:,2),"movmean",window);
y(:,3) = smoothdata(y(:,3),"movmean",window);
%%
k = (k-min(k))./(max(k)-min(k));
k(k<0.5) = 0;
k(k~=0) = 1;

k_R = k(:,1);
k_R(x(:,1)>64 & y(:,1)>100) = 0;
k_R(x(:,1)<64 & y(:,1)>100) = 0;
k_G = k(:,2);
k_G(x(:,2)>64 & y(:,2)>100) = 0;
k_G(x(:,2)<64 & y(:,2)>100) = 0;
k_B = k(:,3);
k_B(x(:,3)>64 & y(:,3)>110) = 0;
k_B(x(:,3)<64 & y(:,3)>110) = 0;
k = [k_R,k_G,k_B];

% get length
idx = find(k(:,1) == 1 & k(:,2) == 1 & k(:,3) == 0);
length_RG = mean(((x(idx,1)-x(idx,2)).^2+(y(idx,1)-y(idx,2)).^2).^0.5);
idx = find(k(:,2) == 1 & k(:,3) == 1 & k(:,1) == 0);
length_GB = mean(((x(idx,2)-x(idx,3)).^2+(y(idx,2)-y(idx,3)).^2).^0.5);
idx = find(k(:,3) == 1 & k(:,1) == 1 & k(:,2) == 0);
length_BR = mean(((x(idx,3)-x(idx,1)).^2+(y(idx,3)-y(idx,1)).^2).^0.5);

%
idx = find(k(:,1) == 0 & k(:,2) == 1 & k(:,3) == 1);
for i = 1:length(idx)
    syms xx yy
    eq1 = (xx - x(idx(i),2))^2 + (yy - y(idx(i),2))^2 == length_RG^2;
    eq2 = (xx - x(idx(i),3))^2 + (yy - y(idx(i),3))^2 == length_BR^2;
    sol = solve([eq1, eq2], [xx, yy]);
    % 提取解
    xx_vals = double(sol.xx); % 转为数值型
    yy_vals = double(sol.yy);
    xy_idx = find(yy_vals>100);
    if isempty(xy_idx)
        k(idx(i),2) = 0;
        k(idx(i),3) = 0;
    else
        x(idx(i),1) = xx_vals(xy_idx);
        y(idx(i),1) = yy_vals(xy_idx);
    end
    % [~,xy_idx] = max(yy_vals);
    % [~,xy_idx] = max([max(xx_vals(1),yy_vals(1));max(xx_vals(2),yy_vals(2))]);
    
end
idx = find(k(:,2) == 0 & k(:,1) == 1 & k(:,3) == 1);
for i = 1:length(idx)
    syms xx yy
    eq1 = (xx - x(idx(i),1))^2 + (yy - y(idx(i),1))^2 == length_RG^2;
    eq2 = (xx - x(idx(i),3))^2 + (yy - y(idx(i),3))^2 == length_GB^2;
    sol = solve([eq1, eq2], [xx, yy]);
    % 提取解
    xx_vals = double(sol.xx); % 转为数值型
    yy_vals = double(sol.yy);
    xy_idx = find(yy_vals>100);
    if isempty(xy_idx)
        k(idx(i),1) = 0;
        k(idx(i),3) = 0;
    else
        x(idx(i),2) = xx_vals(xy_idx);
        y(idx(i),2) = yy_vals(xy_idx);
    end
    % [~,xy_idx] = max(yy_vals);
    % [~,xy_idx] = max([max(xx_vals(1),yy_vals(1));max(xx_vals(2),yy_vals(2))]);
end
idx = find(k(:,3) == 0 & k(:,1) == 1 & k(:,2) == 1);
for i = 1:length(idx)
    syms xx yy
    eq1 = (xx - x(idx(i),1))^2 + (yy - y(idx(i),1))^2 == length_BR^2;
    eq2 = (xx - x(idx(i),2))^2 + (yy - y(idx(i),2))^2 == length_GB^2;
    sol = solve([eq1, eq2], [xx, yy]);
    % 提取解
    xx_vals = double(sol.xx); % 转为数值型
    yy_vals = double(sol.yy);
    xy_idx = find(yy_vals>100);
    if isempty(xy_idx)
        k(idx(i),1) = 0;
        k(idx(i),2) = 0;
    else
        x(idx(i),3) = xx_vals(xy_idx);
        y(idx(i),3) = yy_vals(xy_idx);
    end
    % [~,xy_idx] = max(yy_vals);
    % [~,xy_idx] = max([max(xx_vals(1),yy_vals(1));max(xx_vals(2),yy_vals(2))]);
end

loop = 6;
fy = repmat(fy,loop,1);
fx = repmat(fx,loop,1);

judgement = ones(length(k),1);
idx = find(sum(k,2) ~= 2);
judgement(idx) = 0;
idx = find(k(:,1) == 0 & k(:,2) == 1 & k(:,3) == 1);
judgement(idx) = 0;

idx = find(judgement == 0);
x(idx,:) = [];
y(idx,:) = [];
im(idx,:) = [];
fy(idx,:) = [];
fx(idx,:) = [];
%% get RGB
R_x = x(:,1);
G_x = x(:,2);
B_x = x(:,3);
R_y = y(:,1);
G_y = y(:,2);
B_y = y(:,3);

%% get X Y theta
x = (G_x+B_x+R_x)/3;
y = (G_y+B_y+R_y)/3;
dx = round(x-N/2);
dy = round(y-M/2);

theta = atan2(G_y-B_y,G_x-B_x);
theta = unwrap(theta);
theta = 180/pi*(theta - theta(1));

r_im = im.*exp(1i*2*pi*(fx.*dx+fy.*dy));

%% recover
U = zeros(M,N,3);
for i = 1:length(r_im)
    p = exp(1i*2*pi*(fx(i)*X+fy(i)*Y));
    p = imrotate(p,theta(i),"bilinear",'crop');
    U = U+p.*cat(3,r_im(i,1),r_im(i,2),r_im(i,3));
end
U = abs(U);
U = 255*(U-min(min(U)))./(max(max(U))-min(min(U)));
U = uint8(U);
% run TVAL3
patterns = zeros([length(r_im),M*N]);
for i = 1:length(r_im)
    p = exp(-1i*2*pi*(fx(i)*X+fy(i)*Y));
    p = imrotate(p,theta(i),"bilinear",'crop');
    patterns(i,:) = reshape(p,[1,M*N]);
end

clear opts
opts.mu = 2^4;
opts.beta = 2^13;
opts.tol = 1E-8;
opts.maxit = 1080;
opts.TVnorm = 1;
opts.nonneg = true;

[U_R, ~] = TVAL3(patterns,r_im(:,1),M,N,opts);
[U_G, ~] = TVAL3(patterns,r_im(:,2),M,N,opts);
[U_B, ~] = TVAL3(patterns,r_im(:,3),M,N,opts);
U_TV = cat(3,U_R,U_G,U_B);
U_TV = 255*(U_TV-min(min(U_TV)))./(max(max(U_TV))-min(min(U_TV)));
U_TV = uint8(U_TV);
%% imshow
figure()

tiledlayout(4,4)

nexttile(1,[2,2])
title(data_name)
hold on
plot(R_x,R_y,'.','Color',[1,0,0])
plot(G_x,G_y,'.','Color',[0,1,0])
plot(B_x,B_y,'.','Color',[0,0,1])
plot(x,y,'Color',[0.31764,0.29411,0.29411])
axis square
hold off

nexttile(3,[1,2])
hold on
plot(x,'LineWidth',1.1)
plot(y,'LineWidth',1.1)
hold off
nexttile(7,[1,2])
plot(theta,'LineWidth',1.1)

nexttile(9,[2,2])
imshow(U,[])

nexttile(11,[2,2])
imshow(U_TV,[])
