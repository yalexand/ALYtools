
function [Phi,dPhi,Ind] = adaex_2exp_FLIMFRET(alpha)

    global t; % picoseconds

    global DT;
        
    global IRF;  
    
    global Tp;
    
    global y;
        
    N = numel(t);

    Ndecays = length(y)/N;

    tau1 = alpha(1);
    tau2 = alpha(2);

    A1 = 1/tau1*(1-exp(-Tp/tau1)).*exp(-t/tau1);
    A2 = 1/tau2*(1-exp(-Tp/tau2)).*exp(-t/tau2);
    %
    A1 = DT*conv([ A1'; A1' ],IRF)';
    A1 = A1(N+1:2*N);
    A2 = DT*conv([ A2'; A2' ],IRF)';
    A2 = A2(N+1:2*N);
            
    Phi = zeros(N*Ndecays,2*Ndecays);
    
    for k=1:Ndecays
        index = k*2-1;
        Phi((k-1)*N+1:k*N,index) = A1;
        Phi((k-1)*N+1:k*N,index+1) = A2;        
    end

    Phi = sparse(Phi);    
    
    dPhi = [];
    Ind = [];

     % derivatives - this is correct but doesn't work for "sparse" Phi by
     % some reason
     %
%      % http://www.wolframalpha.com/input/?i=d%2Fdx(1%2Fx%2F(1-exp(-T%2Fx))*exp(-t%2Fx))
%      dA11 = exp((Tp-t)/tau1).*( (t-tau1).*exp(Tp/tau1) - t + Tp + tau1 )/tau1^3/(exp(Tp/tau1)-1)^2;
%      dA22 = exp((Tp-t)/tau2).*( (t-tau2).*exp(Tp/tau2) - t + Tp + tau2 )/tau2^3/(exp(Tp/tau2)-1)^2;
%      %
%      dA11 = DT*conv([ dA11'; dA11' ],IRF)';
%      dA11 = dA11(N+1:2*N);
%      dA22 = DT*conv([ dA22'; dA22' ],IRF)';
%      dA22 = dA22(N+1:2*N);
%      
%      dPhi = zeros(N*Ndecays,2*Ndecays);
%      %
%      for k=1:Ndecays
%         index = k*2-1;
%         dPhi((k-1)*N+1:k*N,index) = dA11;
%         dPhi((k-1)*N+1:k*N,index+1) = dA22;        
%      end
%     
%      Ind = [ 1:2*Ndecays; repmat([1 2],1,Ndecays) ];
%      
%      dPhi = sparse(dPhi);
    
end




















