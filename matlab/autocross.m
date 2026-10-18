% Auto / cross average spectrum
%
% [Output, mediacomp, frequencies] = autocross (data1, data2, fsamp, SecPoints, Win);
%
% data1, data2: time histories input
% fsamp: sampling frequency
% SecPoints: Number of points in each subrecords that it has been used to divide the time history
% Win: time window used to weight the data

function [autocross_mean,frequencies]=autocross(data1,data2,fsamp,SecPoints,Win)


N=length(data1);
df=fsamp/SecPoints;

if (SecPoints/2)==(floor(SecPoints/2))
    frequencies=0:df:(SecPoints/2*df);
else
    frequencies=0:df:((SecPoints-1)/2)*df;
end

NF=length(frequencies);

num_records=fix(N/SecPoints);

autocross=zeros(NF,num_records);

counter=1;
NumPoint=0;

while NumPoint <= N

    start_p=(counter-1)*SecPoints+1;
    
    finish_p=start_p+(SecPoints-1);
    
    
    sp1=fft(Win.*data1(start_p:finish_p));
    sp1=sp1./SecPoints;
    
    sp2=fft(Win.*data2(start_p:finish_p));
    sp2=sp2./SecPoints;
   
    autocross(:,counter)=conj(sp1(1:NF)).*sp2(1:NF);
    
counter=counter+1;

NumPoint=finish_p+SecPoints;

end

if (SecPoints/2)==(floor(SecPoints/2))
    autocross(2:end-1,:)=autocross(2:end-1,:)*2;
else
    autocross(2:end,:)=autocross(2:end,:)*2;
end


autocross_mean=mean(autocross(1:NF,:),2);