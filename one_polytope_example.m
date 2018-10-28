format long
clc
clear all
close all
P=[-2 1; 2 1; 1 2];
fseqG=Gilbert_poly_seq(P);
fseqN=Proj_OnePolytope_seq(P);

%%
figure
semilogy(fseqG(1:50),'b--o', 'LineWidth',2);
hold on;
semilogy(fseqN(1:50),'c--*','LineWidth',1.5);
ylabel({'$\|x_k\|-\|x^*\|$'},'Interpreter','latex','FontSize',13)
xlabel({'Iterations'},'Interpreter','latex','FontSize',13)
legend('Gilbert','NESMINO')

%% Gilbert's algorithm.
function fseqG=Gilbert_poly_seq(P)
z=(1/2)*(P(2,:)'+P(3,:)');
fseqG=zeros(500,1);

for k=1:500
    [a,~]=supp_poly(P,-z);
    z=proj_segment3(a,z);
    %     fseqG(k)=norm(z-[0;1]);
    fseqG(k)=norm(z)-1;
end
end

%% NESMINO
function fseqN=Proj_OnePolytope_seq(P)
mu=0.1;
kappa=0.5;
n=size(P,2);

L=(1/mu)*(norm(P)^2) + 1/2;
alpha=(sqrt(L)-sqrt(kappa))/(sqrt(L)+sqrt(kappa));

u=zeros(n,1);
v=zeros(n,1);
nabla=P'*projsplx(P*v/mu) +(1/2)*v;
fseqN=zeros(500,1);

for i=1:500
    uold=u;
    u=v-(1/L)*nabla;
    v=u + alpha*(u-uold);
    nabla=P'*projsplx(P*v/mu) +(1/2)*v;
    temp=P'*projsplx(P*u/mu);
    fseqN(i)=abs(norm(temp)-1);
end
end
