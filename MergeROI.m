%Merge_ROIs
counter=1;
for fileNb=sorted_index'
    tempidx=sum(idx_corr(counter:fileNb));
    temp=DF(counter:tempidx,:);
    CorrMatrix=zeros(size(temp,1));
    for x=1:size(temp,1)
        for y=1:size(temp,1)
            if y>x
                CorrCoef=corrcoef(temp(x,:),temp(y,:));
                CorrMatrix(x,y)=CorrCoef(1,2);
            end
        end
    end
    
      
    counter=tempidx;