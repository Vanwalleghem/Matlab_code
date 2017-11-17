function DF = DeltaF(RawFluo)
DF=zeros(size(RawFluo));
    parfor x=1:size(RawFluo,1)
        DF(x,:)=(RawFluo(x,:)-mean(RawFluo(x,1:10)))/mean(RawFluo(x,1:10));
    end
end