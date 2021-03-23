clear;
close all;

load data/p_W_landmarks.txt
load data/K.txt
load data/keypoints.txt

N = length(keypoints);
p = [keypoints(:, [2,1])'; ones(1,N)];
P = [p_W_landmarks'; ones(1,N)];

M=K*[eye(3), zeros(3,1)];
p_proj = M*P;

norm(p(1:2,:)./p(3,:)-p_proj(1:2,:)./p_proj(3,:))