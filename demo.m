%This script does the VRBC-t-TNN on a coronary angiogram sequence

%paths and parameters
name = 'wxz2x2';
filename = ['input/' name '.avi'];
outdir = ['output/' name '_3' '/'];
if ~exist(outdir,'dir')
    mkdir(outdir);
end

%read data
vidobj = VideoReader(filename);
n1=vidobj.Height;
n2=vidobj.Width;
n3=vidobj.NumberOfFrames;
vidobj = VideoReader(filename);
D_ori = zeros(n1,n2,n3);
for k = 1 : n3
    frame = readFrame(vidobj);
    if ndims(frame)==3
        frame = rgb2gray(frame);
    end
    frame = double(frame)/255;
    D_ori(:,:,k) = frame;
end

%Preliminary vessel layer extraction
tic
D_log =-log(max(D_ori,0.001));
M = reshape(D_log,n1*n2,n3);
[L,S] = inexact_alm_rpca(M);
L = reshape(L,n1,n2,n3);
S = reshape(S,n1,n2,n3);
% S = D_log-L;
% S(S<0)=0;
S1 = exp(-S);
toc;


%Vessel segmentation
%Radon-Like features
rlfdata = zeros(n1,n2,n3);
tic
for i = 1:n3
    im_fg = S1(:,:,i);
    im_rlf = rlf(im_fg,0.01);
    rlfdata(:,:,i)=im_rlf;
end
%thresholding segmentation and morphological operation
D_vs = zeros(n1,n2,n3,'logical');
threlevel2 = graythresh(S1);
for i = 1:n3
    im_fg = S1(:,:,i);
    im_rlf = rlfdata(:,:,i);
    im_vs = phansalkar(im_rlf,[16,16],1);
    im_vs = ~fillsmallholes(~im_vs,300);
    im_vs2 = ~imbinarize(im_fg,threlevel2);
    im_vs2(im_vs) = 1;
    im_vs0 = im_vs2;
    while ~isequal(im_vs0,im_vs)
        im_vs0 = im_vs;
        im_vs = imdilate(im_vs,ones(3));
        im_vs(~im_vs2) = 0;
    end
    D_vs(:,:,i) = imdilate(im_vs,ones(5));
    imwrite(im_vs,[outdir,'vs_',num2str(i),'.png']);
end
toc;
%Background completion using t-TNN
Omega = ~D_vs;
Idx = find(Omega);
addpath t-TNN;
normalize              =        max(D_log(:))                     ;
Xn                     =        D_log/normalize                   ;
[n1,n2,n3]             =        size(Xn)                      ;
alpha                  =        1.05                             ;
maxItr                 =        100                          ; % maximum iteration
rho                    =        0.01                          ;
A                      =        diag(sparse(double(Omega(:)))); % sampling operator
b                      =        A * Xn(:)                     ; % available data
bb                     =        reshape(b,[n1,n2,n3]);
tic;
[X] = LtSVD_TC(A ,b,rho,alpha ,[n1,n2,n3],maxItr, Xn(:), false, false);
X                      =        X * normalize                 ;
D_bc                  =        reshape(X,[n1,n2,n3]);
D_bc(D_bc<0)=0;D_bc(D_bc>D_log)=D_log(D_bc>D_log);
L0=D_log;L0(D_vs) = D_bc(D_vs);
L0(L0<0)=0;L0(L0>D_log)=D_log(L0>D_log);
S0 = D_log - L0;
eL0 = exp(-L0);eS0 = exp(-S0);
toc;
%save output
for k = 1 : n3
    im = eL0(:,:,k);
    imwrite(im,[outdir,'L0_',num2str(k),'.png']);
    im = eS0(:,:,k);
    imwrite(im,[outdir,'S0_',num2str(k),'.png']);
end