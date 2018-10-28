function P = create_polytopes(p,m,n)
% Generate $p$ polytope of $m$ verteices in $n$ dimension
P = [];
for i = 1 : p
    P = [P;(2*i)*rand(m,n)];
end
end

