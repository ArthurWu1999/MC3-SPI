clc;clear;

%%
M = 128;
N = 128;
data_name = "imR_BIT_RGB_150_TVhad_128x128_6_f = 0.5_1222-235412";
f = 0.5;
%% order
[X, Y] = meshgrid(1:N, 1:M);
values = X+Y+X.*Y;

[~, I] = sort(values(:));
[row_s, col_s] = ind2sub(size(values), I);

%% Data Processing
data = importdata(strcat(data_name,".txt"));

y = (2*data(1:8:end,:)-data(2:8:end,:)-data(3:8:end,:))...
        +1i*3^0.5*(data(2:8:end,:)-data(3:8:end,:));
x = (2*data(4:8:end,:)-data(5:8:end,:)-data(6:8:end,:))...
        +1i*3^0.5*(data(5:8:end,:)-data(6:8:end,:));
im = data(7:8:end,:)-data(8:8:end,:);
r_im = im-mean(im);
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
window = 3;  % 设置高斯平滑的标准差
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
theta = 180/pi*(theta - 0);

%% get patterns
U = zeros(M,N);
num = length(r_im);
patterns = zeros(num,M*N);
for i = 1:num
    s = zeros(M,N);
    s(row_s(i), col_s(i)) = 1;
    p = ifwht(ifwht(s).').';
    p = circshift(p,[-dy(i),-dx(i)]);
    p = imrotate(p,theta(i),"bicubic",'crop');
    U = U+p.*cat(3,r_im(i,1),r_im(i,2),r_im(i,3));
    patterns(i,:) = reshape(p,[1,M*N]);
end
U = 255*(U-min(min(U)))./(max(max(U))-min(min(U)));
U = uint8(U);
% run TVAL3
clear opts
opts.mu = 2^5;
opts.beta = 2^11;
opts.tol = 1E-5;
opts.maxit = 1080;
opts.TVnorm = 1;
opts.nonneg = true;

[U_R, ~] = TVAL3(patterns,r_im(:,1),M,N,opts);
[U_G, ~] = TVAL3(patterns,r_im(:,2),M,N,opts);
[U_B, ~] = TVAL3(patterns,r_im(:,3),M,N,opts);
U_TV = cat(3,U_R,U_G,U_B);
U_TV = 255*(U_TV-min(min(U_TV)))./(max(max(U_TV))-min(min(U_TV)));
U_TV = uint8(U_TV);
% figure();imshow(U,[]);

%%
U_non = zeros(M,N);
U_xy = zeros(M,N);
patterns = zeros(num,M*N);
for i = 1:num
    s = zeros(M,N);
    s(row_s(i), col_s(i)) = 1;
    p = ifwht(ifwht(s).').';
    U_non = U_non+p.*cat(3,r_im(i,1),r_im(i,2),r_im(i,3));
    p = circshift(p,[-dy(i),-dx(i)]);
    U_xy = U_xy+p.*cat(3,r_im(i,1),r_im(i,2),r_im(i,3));
end
U_xy = 255*(U_xy-min(min(U_xy)))./(max(max(U_xy))-min(min(U_xy)));
U_xy = uint8(U_xy);

U_non = 255*(U_non-min(min(U_non)))./(max(max(U_non))-min(min(U_non)));
U_non = uint8(U_non);


%% imshow
figure()
tiledlayout(3,4)

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

nexttile(9)
imshow(U_non,[])

nexttile(10)
imshow(U_xy,[])

nexttile(11)
imshow(U,[])

nexttile(12)
imshow(U_TV,[])
