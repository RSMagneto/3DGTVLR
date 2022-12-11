function [ind_edge,ind_smooth]=find_block(A,sha)
AM=zeros(size(A,3),3);
for i=1:size(A,3)
    AM(i,1)=A(1,1,i);
    AM(i,2)=A(2,2,i);
    AM(i,3)=A(3,3,i);
end
S=AM(:,2).*AM(:,3);
ind_edge=find(S>=sha);
ind_smooth=find(S<sha);

