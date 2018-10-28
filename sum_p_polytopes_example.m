% This program generates $p$ polytopes randomly, each of them is the
% convex hull of $m$ points in $n$ dimension space and then copute
% the projection of the origin to the sum of these p polytopes.

% We compare NESMINO and Gilbert 
%% Initialization
clc;
clear all;
close all;

p=2;        % Number of polytopes,      2<p<=10
m=100;       % Number of vertex of each polytope
n=100;      % Space dimension

prob=10;    % Number of problems to be solved for taking average
rateS=zeros(prob,1001);
rateG=zeros(prob,1001);


for pr=1:prob
    P=create_polytopes(p,m,n);

    mu=10;                                 % for n=500 take mu=100  
    kappa=0.5;                             % strong convexity parameter
    sigma=0.5;                             % for n=100 take sigma=0.5, for n=500 take sigma=0.1,   
    
    %% Use Gilbert's algorithm to find an aprroximate solution for test
    Iter_G=2000;
    sol=Gilbert(P,Iter_G,p);
    %% NESMINO
    u=zeros(n,1);
    v=zeros(n,1);
    normP=norm(P)^2;
    ct=1;
    while ct<201
        L=(1/mu)*normP + 1/2;
        alpha=(sqrt(L)-sqrt(kappa))/(sqrt(L)+sqrt(kappa));
        for k=1:1e4
           nabla=grad_polytope(P,v,mu,p);
            if norm(nabla)> 1e-3
                u_old=u;
                u=v-(1/L)*nabla;
                v=u + alpha*(u-u_old);
                
                ps=primal_seq(P,u,mu,p);
                rateS(pr,ct) = abs(norm(ps) -norm(sol))/norm(sol);
                ct=ct+1;
            else
                break
            end
        end
        mu=sigma*mu;
    end
    
    %% Gilbert's algorithm
    u=rand(n,1);
    z=support_polytope(P,u,p);
    for k=1:201
        a=support_polytope(P,-z,p);
        z=projsegment(a,z);
        rateG(pr,k)=abs(norm(z) - norm(sol))/norm(sol);
    end
    
end
%% Plot
meanrateS=mean(rateS,1);
meanrateG=mean(rateG,1);
figure
semilogy(meanrateG(2:201),'--*r');
hold on;
semilogy(meanrateS(2:201),'--ob');
ylabel({'$\frac{\|x_{k}\| - \|x^*\|}{\|x^*\|}$'},'Interpreter','latex','FontSize',18)
xlabel({'Iterations'},'Interpreter','latex','FontSize',13)
legend('Gilbert','NESMINO')
title(sprintf('$p=%d,m=%d,n=%d$',p,m,n),'Interpreter','latex','FontSize',13)
