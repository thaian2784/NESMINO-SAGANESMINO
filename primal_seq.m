function ps=primal_seq(P,v,mu,p)
ps=P'*reshape(projsplx_matrix(reshape((1/mu)*P*v,[],p)),[],1);
end




