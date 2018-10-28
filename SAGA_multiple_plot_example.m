% SAGA-NESMINO and Gilbert for large scale problems

%% Initialization
format long;
clc;
clear all;
close all;
%%
p=500;           % Number of polytopes, p>>2
m=100;           % Number of vertex of each polytope
n=100;           % Space dimension
fprintf('Generating p=%d polytopes, each has m=%d vertices, in an Euclidean space of n=%d dimension, please wait ...\n\n', p,m,n);

% P=create_polytopes(p,m,n);
load('matrixP_SAG_p500m100n100.mat')       % set p=500, m=100, n=100 before load this
%load('matrixPforSAG_p1000m100n100.mat')   % ok
%%
fprintf('Finding nearly optimal solution x^*, please wait ...\n')
max_iterG=200;
tic
sol=Gilbert(P,max_iterG,p);                          % Referenced solution
time_fG=toc;
%fprintf('Iteration: %d,     Time: %f4, Norm of sol:   %d\n\n', max_iterG, time_fG, norm(sol))
%% Gilbert
fprintf('Performing Gilbert algorithm ...\n\n');
u=[zeros(n-1,1); -1];
z=support_polytope(P,-u,p);
rateG=zeros(1,100);
rateG(1)=(norm(z)-norm(sol))/norm(sol);
for k=2:100
    a=support_polytope(P,-z,p);
    z=projsegment(a,z);
    rateG(k)=((norm(z)-norm(sol))/norm(sol));
%     fprintf('rateG %d: %d\n', k, rateG(k))
end

prob=20;                                      % Number of SAGA plots
fprintf('We are running SAGA-NESMINO in %d times, ... please wait\n', prob)
rateSAGA=zeros(prob,100);
for nprob=1:prob
    fprintf('Performing SAG-NESMINO no. %d  ... \n', nprob)
    %% SAG
    mu=0.01;
    
    x=zeros(n,1);
    d=zeros(n,1);
    y=zeros(n,p);                                          % save p derivatives
    for i=1:p
        Ai=P((i-1)*m+1:i*m,:);
        temp=Ai'*projsplx(Ai*x/mu)+(1/(2*p))*x;
        y(:,i)=temp;
        d=d+temp;
    end
    temp=zeros(p,1);
    for i=1:p
        temp(i)=(norm(P((i-1)*m+1:i*m,:)))^2;
    end
    maxnormA=max(temp);
    L=(1/mu)*maxnormA+(1/(2*p));

    alpha=2500/(16*L);
%   fprintf('\nPerforming SAG-NESMINO with L = %d, alpha=%d\n', L, alpha)
    for q=1:p*100
        temp=randperm(p);
        i=temp(1);
        Ai=P((i-1)*m+1:i*m,:);
        nabla=Ai'*projsplx(Ai*x/mu)+(1/(2*p))*x;
        x=(1-alpha*(1/(2*p)))*x-alpha*(nabla-y(:,i)+(1/p)*d);
        % x=x-alpha*(nabla-y(:,i)+(1/p)*d);
        d=d-y(:,i)+nabla;
        y(:,i)=nabla;
        if mod(q,p)==0
            z=primal_seq(P,x,mu,p);
            t=q/p;
            rateSAGA(nprob,t)=((norm(z)-norm(sol))/norm(sol));
            %fprintf('rateSAGA %d: %d\n', t, rateSAGA(t))
        end
    end
    
end
%% Plot output
fprintf('\nSemilogy plotting ... please wait\n')
figure(1); clf;
semilogy(rateG(1:100),'.r', 'LineWidth',2);
hold on;
for i=1:prob
    semilogy(rateSAGA(i,1:100),'--b','LineWidth',2);
end
ylabel({'$\frac{\|x_{k}\| - \|x^*\|}{\|x^*\|}$'},'Interpreter','latex','FontSize',18)
xlabel({'Iterations'},'Interpreter','latex','FontSize',13)
legend('Gilbert','SAGA-NESMINO')
title(sprintf('$p=%d,m=%d,n=%d$',p,m,n),'Interpreter','latex','FontSize',13)