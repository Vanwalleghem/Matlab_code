function NbCorr = Nb_Corr(RawFluo,GCaMP,timepoints)
    NbCorr=zeros(length(RawFluo),1);
    parfor TraceNb=1:length(RawFluo)
        Trace=RawFluo(TraceNb,:);
        for time=timepoints'
            for lag=-2:5
                time_lag=time+lag;
                temp=Trace(time_lag:time_lag+length(GCaMP)-1);
                temp=corrcoef(temp,GCaMP);
                if temp(1,2)>0.75
                    NbCorr(TraceNb)=NbCorr(TraceNb)+1;
                    break
                end
            end
        end
    end
            
    