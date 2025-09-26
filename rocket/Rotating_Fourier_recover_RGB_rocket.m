clc;clear;

%%
M = 128;
N = 128;
data_name = "imR_rocket_RGB_2000_Fourier_128x128_6_f =0.5_1222-233421";
f = 0.5;

sampling_rate_max = 0.5;
sampling_rate_set = 0.2;
%% order
y_center = floor((M-1)/2)+1;
x_center = floor((N-1)/2)+1;

[X, Y] = meshgrid(0:N-1, 0:M-1);
values = (X-x_center).^2+(Y-y_center).^2;
values = ifftshift(values);
[~, I] = sort(values(:));
[row_s, col_s] = ind2sub(size(values), I);

fy = (row_s-1)/M;
fx = (col_s-1)/N;

%% Data Processing
data = importdata(strcat(data_name,".txt"));
num = length(data);
num = round(num*sampling_rate_set/sampling_rate_max);
num = num-mod(num,9);
data = data(1:num,:);

y = (2*data(1:9:end,:)-data(2:9:end,:)-data(3:9:end,:))...
        +1i*3^0.5*(data(2:9:end,:)-data(3:9:end,:));
x = (2*data(4:9:end,:)-data(5:9:end,:)-data(6:9:end,:))...
        +1i*3^0.5*(data(5:9:end,:)-data(6:9:end,:));
im = (2*data(7:9:end,:)-data(8:9:end,:)-data(9:9:end,:))...
        +1i*3^0.5*(data(8:9:end,:)-data(9:9:end,:));

%% locate
y = -angle(y);
y(y<0) = y(y<0)+2*pi;
y = y/2/pi*M/f;

x = -angle(x);
x(x<0) = x(x<0)+2*pi;
x = x/2/pi*N/f;

%% 高斯平滑
% sigma = 1;  % 设置高斯平滑的标准差
% x1 = imgaussfilt(x1, sigma);
% x2 = imgaussfilt(x2, sigma);
% y1 = imgaussfilt(y1, sigma);
% y2 = imgaussfilt(y2, sigma);
window = 8;  % 设置高斯平滑的标准差
R_x = smoothdata(x(:,1),"movmean",window);
G_x = smoothdata(x(:,2),"movmean",window);
B_x = smoothdata(x(:,3),"movmean",window);
R_y = smoothdata(y(:,1),"movmean",window);
G_y = smoothdata(y(:,2),"movmean",window);
B_y = smoothdata(y(:,3),"movmean",window);
% figure();
% plot(x1,y1)
% hold on
% plot(x2,y2)
% hold on
% plot(x,y)
% hold off

%% get X Y theta
x = (G_x+B_x+R_x)/3;
y = (G_y+B_y+R_y)/3;
dx = round(x-N/2);
dy = round(y-M/2);

theta = atan2(G_y-B_y,G_x-B_x);
theta = unwrap(theta);
theta = 180/pi*(theta-pi/2);

%%
num = length(im);

r_im = im.*exp(1i*2*pi*(fx(1:num).*dx+fy(1:num).*dy));
%% recover
U = zeros(M,N,3);
for i = 1:num
    p = exp(1i*2*pi*(fx(i)*X+fy(i)*Y));
    p = imrotate(p,theta(i),"bilinear",'crop');
    U = U+p.*cat(3,r_im(i,1),r_im(i,2),r_im(i,3));
end
U = abs(U);
U = 255*(U-min(min(U)))./(max(max(U))-min(min(U)));
U = (U-128)*1.15+128;
U(U>255) = 255;
U(U<0) = 0;
U = uint8(U);

% run TVAL3
patterns = zeros([num,M*N]);
for i = 1:num
    p = exp(-1i*2*pi*(fx(i)*X+fy(i)*Y));
    p = imrotate(p,theta(i),"bilinear",'crop');
    patterns(i,:) = reshape(p,[1,M*N]);
end

clear opts
opts.mu = 2^6;
opts.beta = 2^8;
opts.tol = 1E-18;
opts.maxit = 1080;
opts.TVnorm = 2;
opts.nonneg = true;

[U_R, ~] = TVAL3(patterns,r_im(:,1),M,N,opts);
[U_G, ~] = TVAL3(patterns,r_im(:,2),M,N,opts);
[U_B, ~] = TVAL3(patterns,r_im(:,3),M,N,opts);
U_TV = cat(3,U_R,U_G,U_B);
U_TV = 255*(U_TV-min(min(U_TV)))./(max(max(U_TV))-min(min(U_TV)));
U_TV = (U_TV-128)*1.15+128;
U_TV(U_TV>255) = 255;
U_TV(U_TV<0) = 0;
U_TV = uint8(U_TV);
% figure();imshow(U_TV,[]);
%% imshow
figure()
tiledlayout(4,4)

nexttile(1,[2,2])
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
