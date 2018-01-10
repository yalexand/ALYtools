function save_decay_to_FLIMfit_5032_compatible_txt_file(fullfilename,t,v)

 fid = fopen(fullfilename,'w');
 
 fprintf(fid,'%s  %f\r\n','Well',490);
 fprintf(fid,'%s  %f\r\n','WellIndex',0);
 fprintf(fid,'%s  %f\r\n','Wavelength',490);
 fprintf(fid,'%s  %f\r\n','Polarization',2);
 fprintf(fid,'%s  %f\r\n','AcquisitionStatus',1);
 fprintf(fid,'%s  %f\r\n','PowerLevel',4608);
 fprintf(fid,'%s  %f\r\n','AcquisitionTime',13.109596);

    % define maximum length of string
    ML = 0;
    for m=1:numel(t)    
        vt = t(m)/1000;
        vv = v(m);    
        L1 = numel(num2str(vt));
        L2 = numel(num2str(vv));
        ML = max(max(L1,L2),ML);                        
    end

    S = [];
    for k=1:ML
        S = strcat(S,'0');
    end

    for m=1:numel(t)

        vt = t(m)/1000;
        vv = v(m);

        ST = S;
        st = num2str(vt);
        ST(1:length(st))=st;
        %
        SV = S;
        sv = num2str(vv);
        SV(1:length(sv))=sv;

        fprintf(fid,'%s\t%s\r\n',ST,SV);

    end

    fclose(fid);

end

