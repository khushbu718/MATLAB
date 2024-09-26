%function funpso1()
function TC = funpso1 (X)
M=X(1);
T=X(2);
%M=2.5319;
%T=1.1322;
K = 20;    % system capacity				
F = 18;    % durining F-policy				
E = 1;      % epsilon				
v = 0.05;   % gamma using M_n				
b = 1;      % feedback probability				
AA = 0.5;   % aplha				
A=0.05;     % L*B*exp(-A*n); using expo Balking				
B=0.9;      % L*B*exp(-A*n); using expo Balking				
%T=1;        % theta(vacation rate)				
%M=3;        % mu(service rate)				
L=1;        % lambda (arrival rate)				
% CH=30;CO=20;CV=10;CM=10; CT=10;CF=15;		
CH=30;CO=20;CV=20;CM=10;CT=5;CF=20;
mu_case='renege';
lambda_case='expo';
%lambda_case='const';

%% Probabilities
%% state i = 0
P0 = zeros(K,1);
for t = 1:K-1    % for loop calcultes the value of P(0,K)
    prod1 = 1;
    for h = 1:t
        L1 = 1-Gsn(Ln(h));
        prod1 = prod1 * L1;
    end
    P0(t) = (Ln(0) / Ln(t)) * prod1;
end

L2 = Ln(K-1);
G1 = Gsn(L2);
L3 = (1 / G1) - 1;
prod2 = 1;
for h = 1:K-1
    L4 = 1-Gsn(Ln(h));
    prod2 = prod2 * L4;
end
L5 = Ln(0) / Ln(K-1);
P0(K) = L5 * L3 * prod2;


%% State i = 1
P1 = zeros(K+1,1);
P1(1) = (Ln(0)) / AA;
P1(2) = (((Ln(0)) / AA) + 1) * ((Ln(0)) / Mn(1));
%P_1,n for n = 2, 3, ..., F-1
%P11 = zeros(1, F-2);

for t = 2:F-1
    P1(t+1) = Cn(t) * (Ln(0) / Mn(t));
end

%P_1,F
for n=F:F
    prod11=1;
    for h=1:K-1
        L66=1-Gsn(Ln(h));
        prod11=prod11*L66;
    end
    Nume11 = (((Ln(K-1)) / (Mn(K-1))) * (Cn(K-1))) + prod11;
    deno11 = (E / (Mn(K)))*(1 + (((Ln(K-1))/(Mn(K-1))) * Ps(K-1)));
    Nume12=(E/Mn(K))*Nume11;
    deno12=1+deno11;
    P1(n+1) = (Cn(F)-(Nume12/deno12))*((Ln(0))/(Mn(F)));
    
end

%P_1,n n=F,F+1...K-1
for n = F+1:K-1
    %P_1,F
    prod11=1;
    for h = 1:K-1
        L66 = 1 - Gsn(Ln(h));
        prod11 = prod11 * L66;
    end
    Nume11 = ((Ln(K-1) / Mn(K-1)) * Cn(K-1)) + prod11;
    deno11 = (E / (Mn(K)))*(1 + (((Ln(K-1))/(Mn(K-1))) * Ps(K-1)));
    Nume12=(E*Ps(n)/Mn(K))*Nume11;
    deno12=1+deno11;
    P1(n+1) = (Cn(n)-(Nume12/deno12))*((Ln(0))/(Mn(n)));
end

for n = K:K
    prod3 = 1;
    for h = 1:K-1
        L6 = 1-Gsn(Ln(h));
        prod3 = prod3 * L6;
    end
    Nume1 = Ln(0) / Mn(K);
    Nume2 = (Ln(K-1) / Mn(K-1)) * Cn(K-1) + prod3;
    Deno1 = 1 + (((Ln(K-1))/(Mn(K-1))) * Ps(K-1)) ;
    Deno2 = E / Mn(K);
    P1(n+1) = (Nume1 * Nume2) / (1 + (Deno1 * Deno2));
    %P1(n+1) = (Cn(n) - left1) * ((Ln(0)) / Mn(K));
end

%     %P_1,K
%     prod4 = 1;
%     for t = 1:n-1
%         prod4 = prod4 * (1 - Gsn(Ln(t)));
%     end
%     Nume3 = (((Ln(K-1)) / Mn(K-1)) * Cn(K-1)) + prod4;
%     Deno3 = 1 + (((Ln(K-1))/(Mn(K-1))) * (Ps(K-1)));
%     Deno4 = E / Mn(K);
%     left2 = Nume3 / (1 + (Deno3 * Deno4));
%     P1(K+1) = left2 * ((Ln(0)) / Mn(K));

% P_2,n n=F...K
P2 = zeros(K+1-F,1);
for t = 1:K+1-F
    P2(t) = (E / M) * P1(K+1);
end

%% Probabilities
p0 = 1/(sum(P0)+sum(P1)+sum(P2)+1);
Q0 = P0*p0;
Q1 = P1*p0;
Q2 = P2*p0;

%% Performance Measures

PI = sum(Q0)+p0;

PSB = sum(Q1)+sum(Q2);

sum10 = 0;
for t=1:K
    sum10 = sum10 + t*Q0(t);
end
sum11 = 0;
for t=1:K
    sum11 = sum11 + t*Q1(t+1);
end
sum12=0;
for t=F:K
    sum12 = sum12 + (t)*Q2(t-F+1);
end
%Nq(kk) = sum11+sum12;
%Nr(kk) = sum10;
Neq = sum11+sum10+sum12;
% L_eff(kk)=Ln(n)*(sum(Q1)+sum(Q0)-Q1(K+1)-Q0(K));
%Ws(kk)=Neq(kk)/(0.3*L_eff(kk));
%         % system throughput(TP)
%         sum20=0;
%         for t=1:K
%             sum20 = sum20 + Mn(t)*Q1(t+1);
%         end
%         sum21=0;
%         for t=F:K
%             sum21 = sum21 + M*Q2(t-F+1);
%         end
%         TP(kk)=sum20+sum21;

% avarage balking (B)
sum30=0;
for t=1:K
    sum30=sum30+(1-exp(-A*t))*(Q0(t)+Q1(t+1));
end
Bb = L*sum30;

% avarage renege(R)
sum40=0;
for t=1:K
    sum40=sum40+(t*v*Q1(t+1));
end
R=sum40;
CL=R+Bb;
%  TC(kk)=CH*Neq(kk)+CO*PSB(kk)+CV*PI(kk)+M*CM; %cost function
TC=CH*Neq+CO*PSB+CV*PI+M*CM+ T*CT+M*CF; %cost function


%  H=[Mi, Neq ,CL];
% assignin('base','PI',PI);
% assignin('base','Neq',Neq);
% assignin('base','PSB',PSB);
% assignin('base','Q0',Q0);
% assignin('base','Q1',Q1);
% assignin('base','Q2',Q2);
% assignin('base','P1',P1);
% %     assignin('base','L_eff',L_eff);
% %     assignin('base','Ws',Ws);
% assignin('base','TC',TC);
% assignin('base','Bb',Bb);
% assignin('base','R',R);
% assignin('base','CL',CL);
% assignin('base','Z',Z);
% assignin('base','Result',Result);
 assignin('base','TC',TC);

%% Functions
function mn = Mn(n)
% this function give the mu value
switch mu_case
    case 'renege'
        if n >= 1 && n<=K-1
            mn = M + (n - 1) * v;
        elseif n>=K
            mn = M;
        end
end
end

function ln = Ln(n)
% this function finds the value of lambda for 'const' and'
% 'exponential' balking
if n>=0 && n<=K-1
    switch lambda_case
        case 'expo'
            ln = L*exp(-A*n);
        case 'const'
            ln = L*B;
    end
elseif n>=K
    ln = 0;
end
end

function gs = Gsn(t)
% this function will calculate the value of b*(lambda_i)
% t as a input parameters is lambda value

gs = T/(t+T);
% gs = ((3*T)/(t+(3*T)))^3;
%gs = exp(-t/T);
end

%% S_n notation
function c = Cn(n)
% this function calculates the values of f1 eqation
prod1=1;
for ii=1:n-1
    prod1 = prod1 * ((Ln(ii)) / (Mn(ii)));
end
s1 = (((Ln(0)) / AA) + 1) * prod1;
sum1=0;
for m = 1:n-2
    prod2=1;
    for j=1:m
        L7 = Ln(j);
        prod2 = prod2 * (1 - Gsn(L7));
    end
    prod3=1;
    for p = m:n-2
        prod3 = prod3 * ((Ln(p+1)) / (Mn(p+1)));
    end
    sum1=sum1+prod2*prod3;
end
s2=sum1;
prod3=1;
for j=1:n-1
    L7 = Ln(j);
    prod3 = prod3 * (1 - Gsn(L7));
end
s3 = prod3;
c = s1 + s2 + s3;
end
function r = Ps(n)
% this fuction calculates f2 equation
sum2=0;
for q = 0:n-1-F
    prod5=1;
    for tt = q:n-1-F
        prod5 = prod5 * (((Ln(F+tt)) / Mn(F+tt)));
    end
    sum2 = sum2 + prod5;
end
r = sum2 + 1;
end

end





















