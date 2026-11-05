%Matlab script for design of scattering-based dynamic interconnection controller 
%System 1: 
m1=2; %m1 is number of inputs
p1=2; %p1 is the number of the outputs
%W1 is supply rate matrix of the size (m1+p1)*(m1+p1), the first m1
%components are inputs, the last p1 components are outputs 
W1=0.5*[zeros(2) eye(2); eye(2) zeros(2)];
%System 2:
m2=1; %m2 is number of inputs
p2=2; %p2 is the number of the outputs
%W2 is supply rate matrix of the size (m2+p2)*(m2+p2), the first m2
%components are inputs, the last p2 components are outputs 
W2=[0 0.5 0; 0.5 -2 0.5; 0 0.5 -1];
%Liveness condition for system 1:
v1=eig(W1);
n1minus = sum(sign(v1) == -1);
if n1minus==p1
    display('Liveness condition is satisfied for System 1')
else 
    display('Warning: Liveness condition is NOT satisfied for System 1')
end
%Liveness condition for system 2:
v2=eig(W2);
n2minus = sum(sign(v2) == -1);
if n2minus==p2
    display('Liveness condition is satisfied for System 2')
else 
    display('Warning: Liveness condition is NOT satisfied for System 2')
end
% Eigenvalues of W1 in descending order
[G1Prelim, L1Prelim]=eig(W1);
l1=diag(L1Prelim);
[sortedl1, idx1] = sort(l1, 'descend'); 
%the corresponding permutation matrix:
PP1 = eye(m1+p1); 
PP1 = PP1(idx1, :); 
%Calculation of G1 and Lambda1
Lambda1=PP1'*L1Prelim*PP1;
G1=G1Prelim*PP1; 
%
% Eigenvalues of W2 in descending order
%
[G2Prelim, L2Prelim]=eig(W2);
l2=diag(L2Prelim);
[sortedl2, idx2] = sort(l2, 'descend'); 
%the corresponding permutation matrix:
PP2 = eye(m2+p2); 
PP2 = PP2(idx2, :); 
%Calculation of G2 and Lambda2
Lambda2=PP2'*L2Prelim*PP2;
G2=G2Prelim*PP2; 
%Realizability Condition for System 1:
G122=G1([(m1+1):(m1+p1)],[(m1+1):(m1+p1)]);
if rank(G122)==p1
    display('Realizability condition is satisfied for System 1')
    else 
    display('Warning: Realizability condition is NOT satisfied for System 1')
end
%Realizability Condition for System 2:
G222=G2([(m2+1):(m2+p2)],[(m2+1):(m2+p2)]); 
if rank(G222)==p2
    display('Realizability condition is satisfied for System 2')
    else 
    display('Warning: Realizability condition is NOT satisfied for System 2')
end
%Assume permutation matrices P1, P2 in ``Design of S1, S2'' are all unit
%matrices
%Local scattering transformations
S1=G1; 
S2=G2; 
%Design of S1#:
S111=S1([1:m1],[1:m1]);
S112=S1([1:m1],[(m1+1):(m1+p1)]);
S121=S1([(m1+1):(m1+p1)],[1:m1]);
S122=S1([(m1+1):(m1+p1)],[(m1+1):(m1+p1)]);
S122inv=inv(S122);
S1Sharp=[(S111-S112*S122inv*S121) S112*S122inv;
   -S122inv*S121 S122inv];
%Design of S2#:
S211=S2([1:m2],[1:m2]);
S212=S2([1:m2],[(m2+1):(m2+p2)]);
S221=S2([(m2+1):(m2+p2)],[1:m2]);
S222=S2([(m2+1):(m2+p2)],[(m2+1):(m2+p2)]);
S222inv=inv(S222);
S2Sharp=[(S211-S212*S222inv*S221) S212*S222inv;
   -S222inv*S221 S222inv];
%Lambda+ and Lambda- from ``Design of S1, S2'':
lambda1plus=diag(Lambda1([1:m1],[1:m1]));
lambda2plus=diag(Lambda2([1:m2],[1:m2])); 
lambdaplus=cat(1,lambda1plus,lambda2plus); 
LambdaPlus=diag(lambdaplus);
lambda1minus=-diag(Lambda1((m1+1):(m1+p1),(m1+1):(m1+p1)));
lambda2minus=-diag(Lambda2((m2+1):(m2+p2),(m2+1):(m2+p2))); 
lambdaminus=cat(1,lambda1minus,lambda2minus); 
LambdaMinus=diag(lambdaminus);
%
%Design of dynamic controller D:
%Dimension n:
n=max([m1+m2;p1+p2]);
%choose delta
delta=0.5;
%choose LambdaStarPlus
LambdaStarPlus=diag(ones(p1+p2,1));
%find theta and LambdaU:
a_diag=diag(LambdaStarPlus+delta*eye(p1+p2));
b_diag=lambdaminus;
theta=max(a_diag./b_diag);
LambdaU=theta*LambdaPlus+delta*eye(m1+m2);
%Weighting matrix W
W=ones([m1+m2,n])/n; 
%LambdaStarMinus:
LambdaStarMinus=max(svd(W'*LambdaU*W))*eye(n);
LambdaStar=diag(cat(1,diag(LambdaStarPlus),diag(-LambdaStarMinus)));
%
A0=-eye(n); %A0 is chosen as negative unit matrix
B0=eye(n,p1+p2); 
BoldP=lyap(A0,eye(n));
W0=[zeros(p1+p2) B0'*BoldP; BoldP*B0 -eye(n)];
% Eigenvalues/vectors of W0 in descending order
[G0Prelim, LdPrelim]=eig(W0);
ld=diag(LdPrelim);
[sortedld, idx0] = sort(ld, 'descend'); 
%the corresponding permutation matrix:
PP0 = eye(n+p1+p2); 
PP0 = PP0(idx0, :); 
%Calculation of G1 and Lambda1
Lambdad=PP0'*LdPrelim*PP0;
G0=G0Prelim*PP0; 
%Gain GammaD
GammaD=sqrt(inv(Lambdad)*LambdaStar);
GCal=G0*GammaD;
GCal11=GCal([1:(p1+p2)],[1:(p1+p2)]);
GCal12=GCal([1:(p1+p2)],[(p1+p2+1):(p1+p2+n)]);
GCal21=GCal([(p1+p2+1):(p1+p2+n)],[1:(p1+p2)]);
GCal22=GCal([(p1+p2+1):(p1+p2+n)],[(p1+p2+1):(p1+p2+n)]);
A=A0+B0*GCal12*inv(GCal22);
B=B0*(GCal11-GCal12*inv(GCal22)*GCal21);
C=inv(GCal22);
D=-inv(GCal22)*GCal21; 
%Supply rate
TildeW1=[eye(m1) zeros(m1,p1); zeros(m2,m1) zeros(m2,p1)]*W1;
TildeW2=[zeros(m1,m2) zeros(m1,p2); eye(m2) zeros(m2,p2)]*W2;
Wcross=[TildeW1 TildeW2];
WChi=[eye(m1) zeros(m1,p1); zeros(m2,m1) zeros(m2,p1)]*W1*[eye(m1) zeros(m1,p1); zeros(m2,m1) zeros(m2,p1)]'+[zeros(m1,m2) zeros(m1,p2); eye(m2) zeros(m2,p2)]*W2*[zeros(m1,m2) zeros(m1,p2); eye(m2) zeros(m2,p2)]';
W12=[theta*WChi theta*Wcross; theta*Wcross' -eye(m1+m2+p1+p2)];
%minimum eigenvalue of upper a certain matrix:
threshold=0.5; 
littlelambdastar=min(eig(WChi+theta*(1/delta)*Wcross*Wcross'));
if littlelambdastar<threshold
    epsilonnode=-littlelambdastar+threshold;
else
    epsilonnode=0;
end 
HATW12=[theta*(WChi+epsilonnode*eye(m1+m2)) theta*Wcross; theta*Wcross' -eye(m1+m2+p1+p2)];
Upsilon12=(2*theta/delta)*WChi+(4*theta^2/delta^2)*Wcross*Wcross'; 
%vector zeta:=[w_1' y_1' w_2' y_2']'
%dimensions dim w_1=m1, dim y_1=p1, dim w_2=m2, dim y_2=p2
BARW12=HATW12([1, 2, 3, 6, 7, 9, 10], [1, 2, 3, 6, 7, 9, 10]); 
%Display the following variables: 
%S1Sharp
%S2Sharp
%A
%B
%C
%D
%W12
%HATW12
%Upsilon12
%W
%BARW12
