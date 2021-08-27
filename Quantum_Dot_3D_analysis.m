function [P,Q,S,R]=Quantum_Dot_3D_analysis(a,b,t)
%a is the output of up focus plane from IDL, b is the output of bottom focus plane from IDL, t is the frame number of release time point recognized by eye check of metamorph output
%P is (x,y,z) position listed over time, R contains net displacement, velocity and release time point.
%Q is the matrix of each step travallength and speed
%S is the metrix of instantaneous speed
%R is the detailed analysis of parameters
a(:,6)=[]; %remove trajectory
b(:,6)=[]; %remove trajectory
C=[a,b];
[m,~]=size(C);
C(t+1:m,:)=[];
D=zeros(t,1);i=1;
while i<=t
    D(i,1)=0.1*i;
    i=i+1;
end
E=[D,C];
[m,~]=size(E);
j=1;i=1;F=[];
while j<=m
    if E(j,3)>0.8 && E(j,3)<3.0 && E(j,4)>0.8 && E(j,4)<3.0 && E(j,8)>0.8 && E(j,8)<3.0 && E(j,9)>0.8 && E(j,9)<3.0
        F(i,:)=E(j,:);i=i+1;j=j+1;
    else
        j=j+1;
    end
end %reomve bad point with poor ¦Ò
F(:,9)=[];F(:,8)=[];F(:,4)=[];F(:,3)=[];
[m,~]=size(F);
P=zeros(m,5);l=1;
while l<=m
    P(l,1)=F(l,1);
    P(l,2)=(F(l,2)-F(l,5))/(F(l,2)+F(l,5));x=P(l,2);
    P(l,3)=F(l,6);P(l,4)=F(l,7);
    P(l,5)=841.29434-543.27112*x-588.25263*x^2-738.10678*x^3+4226.11262*x^4-853.05745*x^5-5692.01979*x^6;
    l=l+1;
end %get (x,y,z) position correlate to time
P(:,2)=[]; %remove I0
[m,~]=size(P);
x=P(1:m,2);y=P(1:m,3);z=P(1:m,4);
x0=mean(x);y0=mean(y);z0=mean(z);
Rg=sqrt(mean((x-x0).^2+(y-y0).^2+(z-z0).^2)); %calculation of radius of gyration
%Q is the matrix of each step travallength and speed
Q=zeros(m-1,1);j=1;
while j<m
    Q(j,1)=sqrt((P(j+1,2)-P(j,2))^2+(P(j+1,3)-P(j,3))^2+(P(j+1,4)-P(j,4))^2);
    j=j+1;
end
%Q is the matrix of each step travallength and speed
i=1;j=1;S=[];
while i<m
    if P(i+1,1)-P(i,1)<0.2
        S(j,1)=P(i,1);
        S(j,2)=sqrt((P(i+1,2)-P(i,2))^2+(P(i+1,3)-P(i,3))^2+(P(i+1,4)-P(i,4))^2)/0.1;
        i=i+1;j=j+1;
    else
        i=i+1;
    end
end
R=zeros(1,7);
N=sqrt((P(m,2)-P(1,2))^2+(P(m,3)-P(1,3))^2+(P(m,4)-P(1,4))^2);
V=N/(P(m,1)-P(1,1));
R(1,1)=N;%net desplacement nm
R(1,2)=t*0.1-20;%fusion time second
R(1,3)=V;%velocity nm/second
R(1,4)=sum(Q(:,1))/1000;%travelength ¦Ìm
R(1,5)=mean(S(:,2))/1000;%Instantaneous speed ¦Ìm/s
R(1,6)=R(1,4)/(P(m,1)-P(1,1));%Average speed ¦Ìm/s
R(1,7)=Rg;%radius of gyration nm
end
