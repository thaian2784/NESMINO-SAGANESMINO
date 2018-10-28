function sol=Gilbert(P,max_iter,p)
n=size(P,2);
u=rand(n,1);
z=support_polytope(P,u,p);
for k=1:max_iter
    a=support_polytope(P,-z,p);
    z=projsegment(a,z);
end
sol=z;
end