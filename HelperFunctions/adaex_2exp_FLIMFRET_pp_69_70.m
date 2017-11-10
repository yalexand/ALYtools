
function [Phi,dPhi,Ind] = adaex_2exp_FLIMFRET_pp_69_70(alpha)

    global t; % picoseconds

    global DT;
        
    global IRF;  
    
    global Tp;
    
    global y;
        
    N = numel(t);

    Ndecays = length(y)/N;

    tau1 = alpha(1);
    tau2 = alpha(2);

    A1 = 1/DT*conv_irf_pp_69_70(tau1);
    A2 = 1/DT*conv_irf_pp_69_70(tau2);
    
%%%%%%%%%%%%%%%%%%% to see if it matches straight mode of convolution
%     %
%     AA1 = 1/tau1*(1-exp(-Tp/tau1)).*exp(-t/tau1);
%     AA2 = 1/tau2*(1-exp(-Tp/tau2)).*exp(-t/tau2);
%     %
%     AA1 = DT*conv([ AA1'; AA1' ],IRF)';
%     AA1 = AA1(N+1:2*N);
%     AA2 = DT*conv([ AA2'; AA2' ],IRF)';
%     AA2 = AA2(N+1:2*N); 
%     %
%     figure(22);semilogy(t,A1,'bo-',t,A2,'r.-',t,AA1,'co-',t,AA2,'m.-');grid on;
%%%%%%%%%%%%%%%%%%% to see if it matches straight mode of convolution    
        
    Phi = zeros(N*Ndecays,2*Ndecays);
    
    for k=1:Ndecays
        index = k*2-1;
        Phi((k-1)*N+1:k*N,index) = A1;
        Phi((k-1)*N+1:k*N,index+1) = A2;        
    end

    Phi = sparse(Phi);    
    
    dPhi = [];
    Ind = [];
    
end




















