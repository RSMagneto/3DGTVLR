function BH = Extract_blocks(H,para)

sz = size(H);
block_sz = para.block_sz;
overlap_sz = para.overlap_sz;
block_num = (sz-block_sz)./overlap_sz+1;
nblocks = para.nblocks;

BH = zeros([block_sz, nblocks]);
idx = 0;
for i = 1:block_num(1)
    for j = 1:block_num(2)
        for k = 1:block_num(3)
        ii = 1 + (i - 1)*(block_sz(1) - overlap_sz(1));
        jj = 1 + (j - 1)*(block_sz(2) - overlap_sz(2));
        kk = 1 + (k - 1)*(block_sz(3) - overlap_sz(3));
        idx = idx + 1;
        BH(:, :, :, idx) = ...
            H(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, kk:kk+block_sz(3)-1);
        end
    end
end











































