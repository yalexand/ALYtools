
function [Phi,dPhi,Ind] = adaex_1exp_FLIM_pp_69_70(alpha)

    global t; % picoseconds

    global DT;
        
    global IRF;  
    
    global Tp;
    
    global y;
        
    N = numel(t);

    Ndecays = length(y)/N;

    tau1 = alpha(1);

%     A1 = 1/tau1*(1-exp(-Tp/tau1)).*exp(-t/tau1);
%     %
%     A1 = DT*conv([ A1'; A1' ],IRF)';
%     A1 = A1(N+1:2*N);
            
    A1 = 1/DT*conv_irf_pp_69_70(tau1);

    Phi = zeros(N*Ndecays,Ndecays);
    
    for k=1:Ndecays
        index = k;
        Phi((k-1)*N+1:k*N,index) = A1;
    end

    Phi = sparse(Phi);    
    
    dPhi = [];
    Ind = [];
    
end




















