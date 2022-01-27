%% Bending Square Cosserat Analysis of Experiment

clear,clc
%% Finding the best fit parameters

sidelength = [58.22 44.15 30.24 16.15]';            %Enter specimen cross-section side lengths from largest to smallest in mm
Measured_E = [15.05 15.08 25.61 30.84]';            %Enter moduli from BVS bending in MPa from big specimen to small

%Input asymptotic E in MPa from compression
asympE = 9.64;
%Input Poisson's ratio
v = 0.030;

%% Loop to obtain lb, N, and beta/gamma
Omega = Measured_E./asympE;                         %Calculating Omega values
N_index = 0;                                        %Setting index for resulting N values                                                           
boverg = 0.028;                                     %beta/gamma guess, change as needed

for N = 0.23:0.01:0.23                             %N values to search through by increments of 0.01 in this case here.
                                                   %Here, this reads as search N values from 0.2 ro 0.23 in increments 
                                                   %of 0.01.
    N_index = N_index + 1;
    
    %Fitting function for Cosserat elastic solids with square cross
    %sections.  Origin of this equation is a handout from Dr.
    %Drugan,"Bending rigidity in a Cosserat square bar" (2016).  In his
    %paper, he denotes side length as 2a.  Here, side length is just x for
    %convenience and consistency with other equations. 
    fitFunc = @(lb,x) 1+24.*((1+(2.*boverg.*v)+v.^2)./(1+v)).*(lb./(x)).^2- ...
        480.*((boverg+v).^2).*((44-38.*v+3.*(N.^2).*(1-v).*(13-9.*v))./...
        ((N.^2).*(1+v).*(22-19.*v))).*((lb./(x)).^4);
    options = fitoptions('Method','NonLinearLeastSquares', ...
        'StartPoint', 1,'Lower',7.4,'Upper',7.4, ...
        'TolX',1e-15);
    
    %Here I am fitting the data and trying to get the characteristic length
    [fitobject,gof,output] = fit(sidelength,Omega,fitFunc,options);
    parameter1 = fitobject.lb;                      %Characteristic length of bending
    parameter2 = boverg;                            %Reassigning beta over gamma for convenience
        
    %% Calculating R squared
    
    %Creating data for best fit parameters for particular case of current N
    fitFunc1 = @(x) 1+24.*((1+(2.*parameter2.*v)+v.^2)./(1+v)).*(parameter1./(x)).^2- ...
        480.*((parameter2+v).^2).*((44-38.*v+3.*(N.^2).*(1-v).*(13-9.*v))./...
        ((N.^2).*(1+v).*(22-19.*v))).*((parameter1./(x)).^4);
    xx = 0.1:0.1:80;
    yy = fitFunc1(xx);
    
    num = numel(Omega);
    sum_obsv = sum(Omega);
    yavg = (1/num)*sum_obsv;

    %Sum of squares of residuals
    Predicted_O = fitFunc1(sidelength);
    diff = @(d) (Omega(d)-Predicted_O(d)).^2;
    summation = diff(1:num);
    SSres = sum(summation);

    %Total sum of squares
    diff1 = @(d) (Omega(d)-yavg).^2;
    summation1 = diff1(1:num);
    SStot = sum(summation1);

    %Calculating R^2
    R2 = 1-SSres/SStot;
    
    %Compiling matrix of results
    R(N_index,1) = N;
    R(N_index,2) = parameter1;
    R(N_index,3) = parameter2;
    R(N_index,4) = R2;
end


%% Making graph

%This finds the largest value in the speficied column (M) and tells me what
%row it lies in (I)
[M,I] = max(R(:,4));        %Row of highest R^2 value 

pointgenerator = @(x) 1+24.*((1+(2.*R(I,3).*v)+v.^2)./(1+v)).*(R(I,2)./(x)).^2- ...
        480.*((R(I,3)+v).^2).*((44-38.*v+3.*(R(I,1).^2).*(1-v).*(13-9.*v))./...
        ((R(I,1).^2).*(1+v).*(22-19.*v))).*((R(I,2)./(x)).^4);

yy1 = pointgenerator(xx);
yy2 = ones(numel(xx));

%Creating figure
figure
hold on
plot(sidelength,Omega,'ko');
plot(xx,yy1,'k-');
plot(xx,yy2,'r--');
xlim([0 70]);
ylim([0 5]);
xlabel('Width (mm)')
ylabel('Relative Stiffness \Omega')
title(['lb = ',num2str(R(I,2)),'mm, N = ',num2str(R(I,1)),', R^2 = ',num2str(R(I,4)),', b/g = ',num2str(R(I,3)),'.'])
hold off