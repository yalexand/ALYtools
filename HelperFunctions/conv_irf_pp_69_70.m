   function D = conv_irf_pp_69_70(tau)
               
        % See thesis page: 69-70
        %
        global IRF;
        global DT;
        global t;
        global Tp;   
        %
        tg=t;
        g=IRF;
        T=Tp;
        dt=DT; 
        %
        rhoi = exp(tg/tau);
        G = g.*rhoi;
        G = cumsum(G);
        G = circshift(G,1);
        G(1) = 0;
        rho = exp(dt / tau);
        
        A = tau.^2/dt * (1-rho)^2/rho;
        B = tau.^2/dt * (dt / tau - 1 + 1/rho);
        
        C = A * G + B * g .* rhoi;
        
        f = 1 / (exp(T/tau)-1);
        
        D = (C + f * C(end))./ rhoi / tau * dt;
                
    end