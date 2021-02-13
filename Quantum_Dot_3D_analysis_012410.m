function [P]=Quantum_Dot_3D_analysis_012410(a,b,t)
%a is the output of up focus plane from IDL, b is the output of bottom focus plane from IDL, t is the release time point recognized by eye check of metamorph output
%P is (x,y,z) position listed over time, R contains net displacement, velocity and release time point.
%Q is the matrix of each step travallength and speed
%S is the metrix of instantaneous speed
%R is the detailed analysis of parameters
%r is the position from each time point to the last time point
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
    if E(j,3)>0.6 && E(j,3)<5.0 && E(j,4)>0.6 && E(j,4)<5.0 && E(j,8)>0.6 && E(j,8)<5.0 && E(j,9)>0.6 && E(j,9)<5.0
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
    P(l,3)=F(l,3);P(l,4)=F(l,4);
    if x>-0.79668 && x<0.84898
        P(l,5)=1240.29+188.44*log(1.64566/(x+0.79668)-1);
        l=l+1;
    else
        P(l,:)=[];
        l=l+1;
    end
end %get (x,y,z) position correlate to time
p=find( P(:,3)==0);P(p,:)=[];
P(:,2)=[]; %remove I0
end