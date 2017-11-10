
function [Phi,dPhi,Ind] = adaex_sk2_FLIMFRET(alpha)

    global t; % picoseconds

    global DT;
        
    global IRF;  
    
    global Tp;
    
    global y;
    
    global E; % efficiency
    
    global E_avr;
    global tau_FRET_avr;
        
    N = numel(t);

    Ndecays = length(y)/N;

    tauD = alpha(1);
    RDA_div_R0 = alpha(2);

    Nb = numel(t);
    Np = numel(E);

    tau = tauD.*(1-E);

    pE = efficiency_distribution(RDA_div_R0,E);

    FRETSUM = zeros(1,Nb);
    
    for p = 1:Np
        FRETSUM = FRETSUM + pE(p)*(1-exp(-DT/tau(p)))*exp(-t/tau(p))/(1-exp(-Tp/tau(p))); % Cliff's method
        %FRETSUM = FRETSUM + pE(p)/tau(p)*exp(-t/tau(p))/(1-exp(-Tp/tau(p)));
    end
    
    A1 = FRETSUM/sum(pE);                
    A2 = (1-exp(-DT/tauD))/(1-exp(-Tp/tauD)).*exp(-t/tauD); % Cliff's 
    %A2 = 1/tauD/(1-exp(-Tp/tauD)).*exp(-t/tauD);
    
    A1 = conv([ A1'; A1' ],IRF)'; % when using Cliff's method, there is no need to multiply convolution by DT
    A1 = A1(N+1:2*N);
    A2 = conv([ A2'; A2' ],IRF)'; % ditto
    A2 = A2(N+1:2*N);

    Phi = sparse(zeros(N*Ndecays,2*Ndecays));
        
    for k=1:Ndecays
        index = k*2-1;
        Phi((k-1)*N+1:k*N,index) = A1;
        Phi((k-1)*N+1:k*N,index+1) = A2;        
    end

    dPhi = [];
    Ind = [];
    %
    E_avr = sum(pE.*E)/sum(pE);
    tau_FRET_avr = sum(tau.*pE)/sum(pE);
    
end
