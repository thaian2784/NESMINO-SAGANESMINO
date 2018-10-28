function sp=support_polytope(P,u,p)
m=size(P,1)/p;
[~,idx] = max(reshape(P*u,[],p));
mat=zeros(m,p);
mat(sub2ind([m,p],idx',(1:p)'))=1;
sp=P'*reshape(mat,[],1);
end
% when m=1 we shold replace line 3 by idx=ones(1,p)