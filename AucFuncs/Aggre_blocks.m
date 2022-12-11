function B = Aggre_blocks(Y,para,sz)

block_sz = para.block_sz;
overlap_sz = para.overlap_sz;
nblocks = para.nblocks;

B = zeros(sz);
B_ind = zeros(sz);

ker1 = block_sz(1);
ker2 = block_sz(2);
ker3 = block_sz(3);
step1 = overlap_sz(1);
step2 = overlap_sz(2);
step3 = overlap_sz(3);

row_blocksx = 1 + (sz(1) - ker1)/step1;
row_blocksy = 1 + (sz(2) - ker2)/step2;
row_blocksz = 1 + (sz(3) - ker3)/step3;
yz_frame = row_blocksy*row_blocksz;
A = ones(ker1,ker2,ker3);
for k = 1:nblocks
    
    temp = Y(:,:,:,k);
    
    x = (ceil(k/yz_frame) - 1) * step1 + 1; % row
    if mod(k,yz_frame)~=0
        yz = mod(k,yz_frame);
        y = (ceil(yz/row_blocksz) - 1) * step2 + 1;
        if mod(yz,row_blocksz)~=0
            t = (mod(yz,row_blocksz) - 1) * step3 + 1;
        else
            t = (row_blocksz - 1) * step3  + 1;
        end
    else
        
        t = (row_blocksz - 1) * step3 + 1;
        y = (row_blocksy - 1) * step2 + 1;
    end
    
    B(x : x+ker1-1, y : y+ker2-1, t : t+ker3-1) = temp + B(x : x+ker1-1, y : y+ker2-1, t : t+ker3-1);
    B_ind(x : x+ker1-1, y : y+ker2-1, t : t+ker3-1) = A + B_ind(x : x+ker1-1, y : y+ker2-1, t : t+ker3-1);
end
B = B./B_ind;




