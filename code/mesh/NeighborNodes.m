function neigh=NeighborNodes(elements,non)
% function neigh=NeighborNodes(elements,non)
% finds the neighbors of the mesh nodes
maxn=40;%estimate of the maximum number of neighbors; just for initial memory allocation
neigh=zeros(non,maxn+1);
for I=1:non,
    a=find(elements(:,1)==I);
    b=find(elements(:,2)==I);
    c=find(elements(:,3)==I);
    d=[a;b;c];
    nod=elements(d,:);
    temp=unique(nod(:));
    neigh(I,1)=length(temp)-1;
    res=temp(find(temp~=I));
    neigh(I,2:neigh(I,1)+1)=res';
    
end
