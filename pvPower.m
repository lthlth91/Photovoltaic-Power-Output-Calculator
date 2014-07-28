function [plotvalue1,betak, betar,val,val2,val3,plotvalue2,suminso,totalinsolation] = pvPower(filename)
%Created by Tianhao Li. 07/21/2013. Georgia Institute of Technology.
% This program serves the purpose of processing the data from tmy 3 for
% atlanta area. It also calculates the insolation for no reorientation, and
% the reoriented insolation in order to take account peak shaving effect.
% in the first part of the program, the code PROCESS ONLY the six hours
% peak electricity demand time, but in the second part the code decides the
% lost of total energy due to peak shaving effects.


[num] = xlsread(filename);     %high level i/o to import the number values in the tmy3 file.
length = 8760;                 % the whole length of the rows of the values in the TMY3 data sheet.
sumvalue2=zeros(1,24);         %initialzed  the values.
sumvalue3=zeros(1,24);
sumvalue4=zeros(1,24);

%this for loop gets the single hour solar irradiance, for example, for 12pm to 1pm, the average 
%irradiance for the whole year is summation
%of the 365 days hour value and sum it
%seperately for every hour in the 24 hours.
%Additionally, the excel doesn't have easy
%function that can sum up the value
%seperately from one whole column.

for time = 1:24
    for i=time:24:length
        
        sumvalue4(time) = sumvalue4(time) + num(i,4);%this is the summation of the fourth column in excel sheet.
    
    end    
end        
   sumvalue4 = sumvalue4./365; % take the daily average.
   irradiance = sumvalue4;   % solar irradiance for the 
   
 %the following part is not related to the changing of the azimuth angle.
for  zenith = 1:24   
    for i=zenith:24:length
        
        sumvalue2(zenith) = sumvalue2(zenith) + num(i,2); %apply same technique to acquire zenith angles.
    end
end        
   sumvalue2 = sumvalue2./365;
   zenith1 = sumvalue2;
   %the end of the unchanged part, since we are not changing the zenith
   %angle at all.
   
   %cosine = cos(zenith1);
%    sinevalue = sin(zenith1);
   zenithpeak = zenith1(12:17);
   cosinepeak = cos((zenithpeak./180)*pi);
   sinepeak   = sin((zenithpeak./180)*pi);
   
for  azimuth = 1:24   
    for i=azimuth:24:length
        
        sumvalue3(azimuth) = sumvalue3(azimuth) + num(i,3);%for azimuth angles.
    
    end    
end  
% this is the end of exel information processing part.

   sumvalue3 = sumvalue3./365; %daily average for azimuth angles.
   azimu = sumvalue3;
   azimupeak = azimu(12:17); %only peak hours concerned in the first part.
   i=1;    %initialize the index indicator.
   for AzimuthChange = 180:1:250 %this for loop ranges from 180 to 250 in order to find the best azimuth angle for shaving effect.
       %the reason of setting 180 to 250 is because 180 is the best for no
       %shaving, around 1 pm, but shaving we concerned about 3pm, which is
       %around 227, 250 goes to the last peak hour, so 250 is large enough
       %to include the optimal choice.
       AzimuthChange1 = AzimuthChange.*ones(1,6);  %change the size of the angle vector to match with 6 hours original solar angles.
       peakvalue(i,:) = irradiance(12:17).*cos(((AzimuthChange1-azimupeak)./180)*pi).*sinepeak; %calculate the socond coefficient for cos(beta).
       %the above expression is simply cos(r-rs)*sin(zenith), which is
       %regared as the coefficient for cos(beta), beta is tilted angle we
       %want to find, and r is the azimuth angle of the panels we want to
       %find.
       i= i+1; % so the final value for index i will be 72 for sure. the size of the peakvalue is 72* 6.
   end
   
   secondfactor = peakvalue ;   % this part changes. because of peakvalue is related to azimuth.
   firstfactor = zeros(71,6);   %initialize first vector, with size match up with 71 azimuth testing and 6 hours peak values.
   
   for index = 1:71
   
       firstfactor(index,:) = irradiance(12:17).* cosinepeak;  %firstfactor is simply cos(beta)*irradiance. 
                                        %we want to includ the irradiance in the term  in order to get the highest insolation.
   
   end
   
    var1 = atan(secondfactor./firstfactor); %using the formula to sum up A*sin(beta)+B*cos(beta), which gives sqrt(A^2+B^2) * cos(beta-alpha);
    alpha = (var1.*180)./pi;
    R = sqrt(secondfactor.^2 + firstfactor.^2);     %trying to find the biggest R with shifting azimuth angle.
   
    beta = alpha(:,1)-1;    %trials for beta starts at one unit value below the lowest alpha value of the six. But here, since we have 71 trials 
                            %to test the azimuth, which in turn gives the
                            %alpha, so our alpha have 71*6 values. and the
                            %first column has the lowest alpha values.
    beta3 = beta;           
    
    finalvalue = zeros(71,45); 

    for i=1:45 % try 45 times due the range situation of the alpha angles.
        beta = beta +1;
        beta5 = (beta*pi)/180;
        beta2=[beta5 beta5 beta5 beta5 beta5 beta5];  %here beta 
        k = cos(beta2 - var1);
        valsum = R.*cos(beta2-var1); %get six hours insolation here.
        k2 = valsum(:,1)+valsum(:,2)+valsum(:,3)+valsum(:,4)+valsum(:,5)+valsum(:,6); %one day insolation amount.
        finalvalue(:,i) = k2;      
    end
plotvalue1 = finalvalue(35,:);
plotvalue2 = finalvalue(:,22);
betak  = beta3(35):1:(beta3(35)+45);
[val2,ind] = max(finalvalue); % 45 maximum values for insolation every column indicating one sample of beta angle.
                     %45 indices for 45 degrees range of beta angles.
[val,ind2] = max(val2); % the largest insolation we can get for this peak shaving effect purpose from 12pm to 6pm.
ind = ind(ind2);   %get the index for the array of azimuth angles and get the best azimuth angle, at the same time using the
            %previous index value we get the best beta angle.
betabest = beta3(ind);  %the best tilt angle for maximizing insolation.
val2 = finalvalue(1,:);%
val3 = finalvalue(1,:);
betar = beta3(1):1:(beta3(1)+45);
[val2,indt] = max(val2); % the val2 value is the one without reorentation for peak shavings.

 

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %This is the start of the second part.  
 %FIRST WE CONSIDER THE TOTAL INSOLATION AND THE CORRESPONDING TILTED
 %ANGLE.
   zenithpeakT = zenith1(6:20); %everything remains same except that this time we take 24 hours into account.
   azimupeakT = azimu(6:20);    %notice that, don't consider 1 to 5 and 21 to 24 is because no insolation or related activities
   irradT = irradiance(6:20);   %happens at all, so value for that are zero or unuseful.
   cosinepeakT = cos((zenithpeakT./180)*pi);
   sinepeakT   = sin((zenithpeakT./180)*pi);
    
   peakvalueT = irradT.*cos(((180-azimupeakT)./180)*pi).*sinepeakT;
   
   secondfactor2 = peakvalueT ;   % this part changes. because of peakvalue is related to azimuth.
   firstfactor2 = irradT.* cosinepeakT;
  
   var2 = atan(secondfactor2./firstfactor2);       % find the alpha angle in the cobmbined experession R*cos(beta-alpha);
   alpha2 = (var2.*180)./pi;
   R1 = sqrt(secondfactor2.^2 + firstfactor2.^2);  %The amplitude for the new cos equation.
    
   betat= min(alpha2)-1;
    beta6 = zeros(1,15);
    valsumx = zeros(177,15);

    for i=1:177 %range changed to the 24 hours diverse alpha angles.
        betat = betat +1;
        betat2 = (betat*pi)/180;
        beta6(1,:)= betat2;  
        k = cos(beta6 - var2);
        valsum2(i) = sum(R1.*cos(beta6-var2));
        valsumx(i,:)= R1.*cos(beta6-var2);
    end
    
       % plotvalues = valsum2(1:60)
        [ totalinsolation indexb] = max(valsum2); % total insolation without changes for shaving effects.
        insolationbyhour = valsumx(indexb,:);
        
        
        
%SECONDLY, SIMPLY CALCULATE THE TOTAL INSOLATION FOR SHAVING EFFECTS, WHICH
%MEANS THE INSOLATION FOR REORIENTED PV SYSTEMS FOR 24 RATHER 6 HOURS.

newbeta = 36.1683; %using the fixed optimal tilted angle calculated from part1.
newbeta = (newbeta*pi)/180; 
newazimuth = 214;   %the optimal azimuth angle for shaving effect.
peakvaluek = irradT.*cos(((newazimuth-azimupeakT)./180)*pi).*sinepeakT; 

secondfactor3 = peakvaluek;    % 
var3 = atan(secondfactor3./firstfactor2);       % find the alpha angle in the cobmbined experession R*cos(beta-alpha);
R3 = sqrt(secondfactor3.^2 + firstfactor2.^2) ; %The amplitude for the new cos equation.

newbeta2 = zeros(1,15);
newbeta2(1,:) = newbeta;
newinsolation = R3.*cos(newbeta2-var3);
suminso = sum(newinsolation);  %total insolation per day for shaving-shifted solar panels.
     
   
   


