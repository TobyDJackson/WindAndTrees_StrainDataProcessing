%MAKE A VECTOR function for this where I take a vector and find those where
%both signals are positive (case1) and the others (case2) and deal with
%them separately in vector operations before combining them back as before.

function [max_strain, angle] = Solve_vector( S1, S2, th1, th2 )
%SG2_bearing must be greater than SG1_bearing! If not, rearrange the channels.
if th1 >= th2 %Switch them around via temporary variables
    temp1=S1; temp2=S2; temp3=th1; temp4=th2;
    S1=temp2; S2=temp1; th1=temp4; th2=temp3;
    clearvars temp1 temp2 temp3 temp4
end

%S1 and S2 are the calibrated and smoothed STRAINS from the two strain
%gauges, -ve values are compression and +ve values are tension.
%The output 'angle' is the angle from S2 to the position of max strain - in
%the direction away from the S1 

separation=abs(th2-th1);
%If they are perpendicular I go straight to trig. If not I correct.
if separation ~=90
    diff=separation-90;
    S2=S2*cosd(abs(diff));
    th2=th2-diff;
end

%Now my data is perpendicular and I can go down two routes
% 1. Simple trigonometry to find the max strain and its angle
% 2. Find the strain on two fixed directions, N & E perhaps. 

% for now I do 1.
%max_strain=sqrt((S1.*abs(S1))+(S2.*abs(S2))); %IS THIS RIGHT?
max_strain=sqrt((S1.^2)+(S2.^2));
angle=th1+atand(S2./S1);
angle(find(S1<0))=angle(find(S1<0))+180;
angle(find(angle>360))=angle(find(angle>360))-360;

end
