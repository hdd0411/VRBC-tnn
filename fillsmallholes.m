function new0=fillsmallholes(bw0,threshold)
% Fill small holes in the binary image bw using the threshold, only holes
% with areas smaller than threshold will be filled
bw = ones(size(bw0)+[2,2]);
bw(2:end-1,2:end-1) = bw0;
filled = imfill(bw, 'holes');

holes = filled & ~bw;

bigholes = bwareaopen(holes, threshold);

smallholes = holes & ~bigholes;

new = bw | smallholes;
new0 = new(2:end-1,2:end-1);