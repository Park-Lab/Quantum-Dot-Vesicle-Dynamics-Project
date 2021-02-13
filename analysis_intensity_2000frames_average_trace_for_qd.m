function [D]=analysis_intensity_2000frames_average_trace_for_qd(Q)
m=length(Q); c=m/2000;i=1;j=1;p=1;% divide the data from 1 Column into a matrix (A) with 200 Row and n Column
A=zeros(2000,c);
while j<=c
    while i<=2000
        A(i,j)=Q(p,1);p=p+1;
        i=i+1;
    end
    j=j+1;i=1;
end
[m,n]=size(A);B=zeros(m,n-1);% remove the background(supposed to be the last Column of matrix A)
v=1;w=1;
while w<=n-1
    while v<=m
        B(v,w)=A(v,w)-A(v,n);v=v+1;
    end
    v=1;w=w+1;
end

k=zeros(1,n-1);
h=1;
while h<=n-1
    k(1,h)=mean(B(1:200,h));% k is the average of baseline (first 20 seconds of each dot); h column, each column (each area) has only 1 averaged baseline from their frist 20 dots
    h=h+1;
end
C=zeros(m,n-1); g=1;h=1;% C is the matrix of data after normalized with baseline
while h<=n-1
    while g<=m
        C(g,h)=B(g,h)/k(1,h);
        g=g+1;
    end
    g=1;h=h+1;
end
[m,n]=size(C);
D=zeros(m,3);g=1;
while g<=m
    D(g,1)=g;g=g+1;
end
g=1;
while g<=m
    D(g,2)=mean(C(g,1:n));g=g+1;
end
g=1;
while g<=m
    D(g,3)=std(C(g,1:n))/sqrt(n);
    g=g+1;
end
end