function [R,A,P]=D3_steering_block_3D(H,para)

sz = size(H);
block_sz = para.block_sz;
overlap_sz = para.overlap_sz;
nblocks = para.nblocks;

R = zeros(3, 3, nblocks);

A = zeros(3, 3, nblocks);

P = zeros(block_sz(1), block_sz(2), block_sz(3), nblocks);

K = fspecial_3D('disk', block_sz(1));

[Hx,Hy,Hz]=gradient(H);

block_num = (sz-block_sz)./overlap_sz+1;
idx = 0;
epslon = 10e-6;
for i = 1 : block_num(1)
    for j = 1 : block_num(2)
        for k = 1 : 1:block_num(3)
            
            ii = 1 + (i - 1)*(block_sz(1) - overlap_sz(1));
            jj = 1 + (j - 1)*(block_sz(2) - overlap_sz(2));
            kk = 1 + (k - 1)*(block_sz(3) - overlap_sz(3));
            
            gx = Hx(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, kk:kk+block_sz(3)-1) .* K;
            gy = Hy(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, kk:kk+block_sz(3)-1) .* K;
            gz = Hz(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, kk:kk+block_sz(3)-1) .* K;
            
            G = [gx(:), gy(:), gz(:)];
            
            [~, s, v] = svd(G, 0);
            % Smallest singular value in s correspondes edge direction in v            
            
            idx = idx + 1;
            
            R(:,:,idx) = v; 
            
            A(:,:,idx) = s / (s(1,1)+epslon);
            
            P(:,:,:,idx) = H(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, kk:kk+block_sz(3)-1);
        end
    end
end






