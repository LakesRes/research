%% Torsion Round Cosserat Analysis of Experiment

clear,clc, close all
%% Finding best fit parameters

%Input the diameter of tested specimens and their tested shear moduli
%These are the 111 structures, G = 3.9 MPa, lt = 6.0 mm
diam = [80.53 44.3 22.21]';             %Enter diameters in mm from largest specimen to smallest
Measured_G = [4.64 5.38 11.1]';         %Enter corresponding shear moduli in MPa, largest specimen to smallest. 

%% For loops to determine best fit G' and lt
%var is asymptotic G in MPa
%lt is mm

count1 = 0;     %Setting counting index for final solutions matrix.
%First cycle through a range of shear moduli.  Specify where to
%start and stop in MPa, to set a single value make start and stop points 
%the same (same for characteristic length, lt, to come). Specify the the 
%interval to search. 
start = 3.9;
stop = 3.9;
interval = 0.1;
for var = start:interval:stop
    count1 = count1 + 1;
    
    %Cycle through a range of characteristic lengths in mm
    count2 = 0;
    startlt = 6.0;
    stoplt = 6.0;
    psi = 1.5;          % Set the value of psi as desired here. 
    for lt = startlt:0.1:stoplt
        count2 = count2 + 1;
        
        %Output is equivalent to relative stiffness 
        Output = Measured_G./var;

        %Function I am fitting the data with.  Origin of this function is
        %Gauthier and Jahsman's 1975 paper, "A Quest for Micropolar Elastic
        %Constants." p leads to N and is defined in this paper. To set a
        %search a range of N values, the equivalent p values must be
        %calculated.
       
        fitFunc = @(p,x) 1+6*((lt/1000)./(x./2000)).^2.*((1-4/3*psi.* ...
            (besseli(1,(p.*x./2000))./((p.*x./2000).*besseli(0,(p.* ...
            x./2000)))))./(1-psi.*(besseli(1,(p.*x./2000))./((p.*x./ ...
            2000).*besseli(0,p.*x./2000)))));
        options = fitoptions('Method','NonLinearLeastSquares',...
            'StartPoint',-10000,'Lower',-10000,'Upper',10000,'TolX',1e-25);

        %Fitting output: p
        %p is a parameter specified in the torsion equation that 
        %is used to calculate torsional elastic constants
        [fitobject,gof,output] = fit(diam,Output,fitFunc,options);
        parameter = fitobject.p;

        %This is a set of test data using the parameters calculated
        %above to determine an R^2 value
        xx = (1:0.1:100)';
        test = @(x) 1+6*((lt/1000)./(x./2000)).^2.*((1-4/3*psi.* ...
            (besseli(1,(parameter.*x./2000))./((parameter.*x./ ...
            2000).*besseli(0,(parameter.*x./2000)))))./(1-psi.* ...
            (besseli(1,(parameter.*x./2000))./((parameter.*x./ ...
            2000).*besseli(0,parameter.*x./2000)))));
        yy=test(xx);
                
        %Calculating R^2
        num = numel(Output);
        sum_obsv = sum(Output);
        yavg = (1/num)*sum_obsv;

        %Sum of squares of residuals
        Predicted_G = test(diam);
        diff = @(d) (Output(d)-Predicted_G(d)).^2;
        summation = diff(1:num);
        SSres = sum(summation);

        %Total sum of squares
        diff1 = @(d) (Output(d)-yavg).^2;
        summation1 = diff1(1:num);
        SStot = sum(summation1);

        %Calculating R^2
        R2(count2,1) = 1-SSres/SStot;
        R2(count2,2) = lt;
        R2(count2,3) = parameter;
    end
    
    %Here I store the highest R^2 value and its corresponding lt
    %from the entire range of tested lt values. This is done for
    %each value of G' (asymptotic G). 
    [maxi idx] = max(R2(:,1));
    Final(count1,1) = maxi;         %max R^2
    Final(count1,2) = R2(idx,2);    %best fit lt
    Final(count1,3) = var;          %G' at which fit occured
    Final(count1,4) = R2(idx,3);    %p required for best fit
    
end

%% Extracting global best fit parameters
[maxR2 indexval] = max(Final(:,1)); %Maximum based on highest R^2
asympG = Final(indexval,3);         %Best asymptotic G
bestfitlt = Final(indexval,2);      %Best characteristic length in torsion
bestp = Final(indexval,4);          %p at best fit

%% Creating a plot to illustrate the best fitting parameters
% range = start:interval:stop;
% hold on
% plot(range,Final(:,1),'k-');
% xlabel('Asymptotic G (MPa)');
% ylabel('R^2');
% hold off
% title(['Best fit when Asymptotic G = ',num2str(asympG),' MPa and lt = ',num2str(bestfitlt),'mm. R^2 = ',num2str(maxR2),'.']);
% hold off

%% Creating the relative stiffness plot
figure
hold on
Omega = Measured_G./asympG;     %Relative stiffness according to best fit
lt = bestfitlt;                 %Redefining lt to be best fit lt
fitFunc1 = @(p,x) 1+6*((lt/1000)./(x./2000)).^2.*((1-4/3*psi.* ...
            (besseli(1,(p.*x./2000))./((p.*x./2000).*besseli(0,(p.* ...
            x./2000)))))./(1-psi.*(besseli(1,(p.*x./2000))./((p.*x./ ...
            2000).*besseli(0,p.*x./2000)))));
yy1 = fitFunc1(bestp,xx);
yy2 = ones(size(xx,1));
plot(diam,Omega,'ko');
plot(xx,yy1,'k-');
plot(xx,yy2,'r--');
xlim([0 90])
ylim([0 5])

%Calculating Cosserat parameters
bplusg = (bestfitlt/1000).^2.*2.*asympG.*10^6;
alpha = (bplusg./1.5)-bplusg;
Kappa = (bestp.^2.*(alpha+bplusg))/2;
N = (Kappa./(2.*asympG*10^6+Kappa)).^(1/2);

%Title and labels have to come after calculation of N
title(['Relative Stiffness vs. Diameter in Torsion. N = ',num2str(N),'. R^2 = ',num2str(maxR2),'.'])
xlabel('Diameter (mm)')
ylabel('Relative Stiffness (\Omega)')
hold off