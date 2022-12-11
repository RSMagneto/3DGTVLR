clear all;
close all;
addpath(genpath(pwd));

% Load data
load('NoiseImage.mat')

% Parameter settings
para.block_sz = [8,8,8];
para.overlap_sz = [4,4,4];
para.alpha = 1; 
para.beta = 0.0042; 
para.lambda = 0.0004;
para.lambda2 = 0.0002;

% Extension for proper size
EH = Pad_HS(NoiseImage,para);
ss = size(EH);
para.nblocks = prod((ss-para.block_sz)./para.overlap_sz+1);

% Direction estimation
[R,A,P_all]=D3_steering_block_3D(EH,para);
Direc.R=R;
Direc.A=A;
Direc.sha=0.8;

% Denoising
[MH,MS,BH] = GTV_HS_Denoising_GloLR_3DGTV_12_lambda12(EH,Direc,para);















