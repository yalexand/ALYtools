
function [Phi,dPhi,Ind] = adaex_2exp_FLIMFRET_fixed_tauD_pp_69_70(alpha)

    global t; % picoseconds

    global DT;
        
    global IRF;  
    
    global Tp;
    
    global y;
    
    global fixed_tauD;
        
    N = numel(t);

    Ndecays = length(y)/N;

    tau1 = fixed_tauD;
    tau2 = alpha(1);

%     A1 = 1/tau1*(1-exp(-Tp/tau1)).*exp(-t/tau1);
%     A2 = 1/tau2*(1-exp(-Tp/tau2)).*exp(-t/tau2);
%     %
%     A1 = DT*conv([ A1'; A1' ],IRF)';
%     A1 = A1(N+1:2*N);
%     A2 = DT*conv([ A2'; A2' ],IRF)';
%     A2 = A2(N+1:2*N);

    A1 = 1/DT*conv_irf_pp_69_70(tau1);
    A2 = 1/DT*conv_irf_pp_69_70(tau2);

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




















