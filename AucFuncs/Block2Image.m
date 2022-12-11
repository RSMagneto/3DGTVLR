function [Y,S,Y_cubeRecon]=Block2Image(Y_all_blocks,Y_all_column_blocks,S_all_column_blocks,para)

s = para.s;
Y=zeros(s);
S=zeros(s);
Y_cubeRecon=zeros(s);
T=zeros(s);

block_sz = para.block_sz;
overlap_sz = para.overlap_sz;
block_num = (s-block_sz)./overlap_sz+1;
T_sub = ones(block_sz);
idx = 0;

for i = 1 : block_num(1)
    for j = 1 : block_num(2)
        for k = 1 : 1:block_num(3)
            
            idx = idx + 1;
            
            ii = 1 + (i - 1)*(block_sz(1) - overlap_sz(1));
            jj = 1 + (j - 1)*(block_sz(2) - overlap_sz(2));
            kk = 1 + (k - 1)*(block_sz(3) - overlap_sz(3));
            
            Y(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, kk:kk+block_sz(3)-1) = ...
                Y_all_column_blocks(:,:,:,idx) + Y(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, kk:kk+block_sz(3)-1);
            
            S(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, kk:kk+block_sz(3)-1) = ...
                S_all_column_blocks(:,:,:,idx) + S(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, kk:kk+block_sz(3)-1);
            
            Y_cubeRecon(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, kk:kk+block_sz(3)-1) = ...
                Y_all_blocks(:,:,:,idx) + Y_cubeRecon(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, kk:kk+block_sz(3)-1);
             
            T(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, kk:kk+block_sz(3)-1) = ...
                T_sub + T(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, kk:kk+block_sz(3)-1);
            
        end
    end
end

Y = Y ./ T;
S = S ./ T;
Y_cubeRecon = Y_cubeRecon ./ T;

