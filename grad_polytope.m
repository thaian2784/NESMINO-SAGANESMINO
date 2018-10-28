function grad=grad_polytope(P,v,mu,p)
grad=P'*reshape(projsplx_matrix(reshape((1/mu)*P*v,[],p)),[],1)+ (1/2)*v;
end




