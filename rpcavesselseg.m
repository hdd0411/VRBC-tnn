clear;clc;
%paths and parameters
name = 'mxl2x';
filename = ['input/' name '.avi'];
outdir = ['output/' name '/'];
if ~exist(outdir,'dir')
    mkdir(outdir);
end

%read data
n1=512;n2=512;n3=80;
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

%log process
D_log =-log(max(D_ori,0.001));

%rpca
M = reshape(D_log,n1*n2,n3);
[L,S] = inexact_alm_rpca(M);
L = reshape(L,n1,n2,n3);
S = D_log-L;
S(S<0)=0;
S1 = exp(-S);

%Radon-Like features
rlfdata = zeros(n1,n2,n3);
for i = 1:n3
    im_fg = S1(:,:,i);
    im_rlf = rlf(im_fg,0.01);
    rlfdata(:,:,i)=im_rlf;
end

%thresholding segmentation and morphological operation
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
    imwrite(im_vs,[outdir,'vs_',num2str(i),'.png']);
end