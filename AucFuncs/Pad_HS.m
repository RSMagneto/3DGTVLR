function EH = Pad_HS(H,para)

[m,n,b] = size(H);
block_sz = para.block_sz;
overlap_sz = para.overlap_sz;

if(mod(m-block_sz(1),overlap_sz(1))~=0)
    disp('X���ص����ܱ�����');
    mm = ceil((m-block_sz(1))/overlap_sz(1))*overlap_sz(1) + block_sz(1);
else
    mm = m;
end

if(mod(n-block_sz(2),overlap_sz(2))~=0)
    disp('Y���ص����ܱ�����');
    nn = ceil((n-block_sz(2))/overlap_sz(2))*overlap_sz(2) + block_sz(2);
else
    nn = n;
end

if(mod(b-block_sz(3),overlap_sz(3))~=0)
    disp('Z���ص����ܱ�����');
    bb = ceil((b-block_sz(3))/overlap_sz(3))*overlap_sz(3) + block_sz(3);
else
    bb = b;
end

EH = padarray (H,[mm-m,nn-n,bb-b],'symmetric','post');














