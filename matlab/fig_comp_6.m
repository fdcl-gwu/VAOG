clear;
close all;

load comp_6

load data/p_W_landmarks.txt;
load data/K.txt;
load data/keypoints.txt;
keypoints = keypoints';




figure(1);plot(1:N,JJ);
hold on;
set(gca,'xscale','log','yscale','log');
legend('$\mathrm{ELGVI}\; (p=3, h=0.00015)$',...
    'Location','SouthWest','interpreter','latex');
ylabel('${f}(R)-{f}(R^*)$','interpreter','latex');
xlabel('$k$','interpreter','latex');
set(gca,'FontSize',16);

print('comp_6a','-depsc');


figure(2);
N = length(keypoints);
p = [keypoints([2,1],:); ones(1,N)];
P = [p_W_landmarks'; ones(1,N)];

M = K*[R0, x0];
p_proj = M*P;

keypoints_proj = p_proj(1:2,:)./p_proj(3,:);
keypoints_proj = keypoints_proj([2 1],:);

img = imread('data/000000.png');
imshow(img);
hold on;
scatter(keypoints(2,:), keypoints(1,:),'r+', 'LineWidth',2);
hold on;
scatter(keypoints_proj(2,:), keypoints_proj(1,:),'b+', 'LineWidth',2);
plot([keypoints(2,:); keypoints_proj(2,:)], [keypoints(1,:); keypoints_proj(1,:)],'b:');

print('comp_6b','-depsc');


figure(3);
M = K*[R, x];
p_proj = M*P;

keypoints_proj = p_proj(1:2,:)./p_proj(3,:);
keypoints_proj = keypoints_proj([2 1],:);

img = imread('data/000000.png');
imshow(img);
hold on;
scatter(keypoints(2,:), keypoints(1,:),'r+', 'LineWidth',2);
hold on;
scatter(keypoints_proj(2,:), keypoints_proj(1,:),'b+', 'LineWidth',2);
plot([keypoints(2,:); keypoints_proj(2,:)], [keypoints(1,:); keypoints_proj(1,:)],'w');


print('comp_6c','-depsc');

!cp comp_6*.eps ../figs