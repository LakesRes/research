%% Torsion Square Cosserat Analysis of Experiment
%N = 1 approximation

clear,clc

%% Finding the best fit parameters

%Entering data and preliminary calculations
sidelength = [58.22 44.15 30.24 16.15]';        %Enter specimen cross-section side lengths from largest to smallest in mm
Measured_G = [4.56 4.18 7.1 9.3]';              %Enter moduli from BVS torsion in MPa from big specimen to small

asymptG = 3.3;                                  %Enter asymptotic G guess in MPa.  Change this to to find R^2 at each guess. 
RelativeStiff = Measured_G./asymptG;            %Calculating Omega values for each G
lb = 3.3;                                       %Characteristic length guess in mm

index = 0;          %Setting index for final matrix organization.

%These function are from Dr. Drugan and Dr. Lakes's 2018 paper, "Torsion of a Cosserat
%elastic bar with square cross section: theory and experiment."  
lbbar = @(lb,x) 2.*lb./(x);                 %Function for characteristic length of bending
lbar = @(lt,x) 2.*lt./(x);                  %Function for characteristic length of torsion
    
%This equation calls the functions for lbbar and lbar to fit for lt and lb
%This equation is for the specific case of N = 1.  Origin of this equation
%is from Dr. Drugan and Dr. Lakes's 2018 paper, "Torsion of a Cosserat
%elastic bar with square cross section: theory and experiment." lbar here
%is the characteristic length of torsion and lbbar is the characteristic
%length of bending.  

N = 1;
    fitFunc = @(lt,lb,x) ((4./21).*((1796 + 126.*(449+2740.*(lbar(lt,x)).^2 ...
        +3960.*(lbar(lt,x)).^4).*(lbar(lt,x)).^2 + 693.*(152+2280.* ...
        (lbar(lt,x)).^2+6615.*(lbar(lt,x)).^4).*(lbbar(lb,x)).^2)./ ...
        (8.*(19+465.*(lbar(lt,x)).^2+990.*(lbar(lt,x)).^4) + 1485.*(6+ ...
        49.*(lbar(lt,x)).^2).*(lbbar(lb,x)).^2)))./(898./399);
 
 %The first term in the square brackets within the option below control lt
 %while the second term controls lb.  These can be changed to force their
 %respective limits as indicated.  
    options = fitoptions('Method','NonLinearLeastSquares','StartPoint', ...
        [1.0*10^(0),3.3*10^(0)],'Lower',[3.8*10^(0),3.3*10^(0)],'Upper', ...
        [1.0*10^(5),3.3*10^(0)],'TolX',1e-10);
    
    [fitobject,gof,output] = fit(sidelength,RelativeStiff,fitFunc,options);
    parameter1 = fitobject.lt;                  %Characteristic length of torsion
    parameter2 = lb;                            %Characteristic length of bending
          
%% Calculating R^2 
    xx = (0.1:0.1:80)';
    yy = fitFunc(parameter1,parameter2,xx);
    
    num = numel(RelativeStiff);
    sum_obsv = sum(RelativeStiff);
    yavg = (1/num)*sum_obsv;

    %Sum of squares of residuals
    Predicted_Omega = fitFunc(parameter1,parameter2,sidelength);
    diff = @(d) (RelativeStiff(d)-Predicted_Omega(d)).^2;
    summation = diff(1:num);
    SSres = sum(summation);

    %Total sum of squares
    diff1 = @(d) (RelativeStiff(d)-yavg).^2;
    summation1 = diff1(1:num);
    SStot = sum(summation1);

    %Calculating R^2
    R2 = 1-SSres/SStot;
    
    %Entering values into table
    R(1,2) = parameter1;            %lt (mm)
    R(1,3) = parameter2;            %lb (mm)
    R(1,4) = R2;                    %R^2

%% Plotting using best fit parameters called from R matrix

[max,idx] = max(R(:,4));
fitFunc2 = @(lt,lb,x) ((4./21).*((1796 + 126.*(449+2740.*(lbar(lt,x)).^2 ...
        +3960.*(lbar(lt,x)).^4).*(lbar(lt,x)).^2 + 693.*(152+2280.* ...
        (lbar(lt,x)).^2+6615.*(lbar(lt,x)).^4).*(lbbar(lb,x)).^2)./ ...
        (8.*(19+465.*(lbar(lt,x)).^2+990.*(lbar(lt,x)).^4)+1485.*(6+ ...
        49.*(lbar(lt,x)).^2).*(lbbar(lb,x)).^2)))./(898./399);
xx1 = (0.1:0.1:100);
yy1 = fitFunc2(R(idx,2),R(idx,3),xx1);
yy2 = ones(1,numel(xx1));

%Creating figure here
figure
hold on
data = plot(sidelength,RelativeStiff,'k.','MarkerSize',30);
fit1 = plot(xx1,yy1,'k-');
fit2 = plot(xx1,yy2,'k--');
xlim([0 70]);
ylim([0 5]);
xlabel('Width (mm)')
ylabel('Relative Stiffness \Omega')
title(['Best fit for torsion. G = ',num2str(asymptG),' MPa, lt = ',num2str(R(idx,2)),'mm, lb = ',num2str(R(idx,3)),'mm, N = ',num2str(N),', R^2 = ',num2str(R(idx,4)),'.']);
%title(['Enter title here during preliminary stages and uncomment']);
hold off
