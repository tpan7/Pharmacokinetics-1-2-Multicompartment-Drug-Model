function [solvedKeK, solvedKeL] = ke_solveMC(parameters, Dose)

flowrates = 60*[4.7 0.750 0.3413 0.0597 0.3917 0.7838 1 0.705]; %L/hr
volumes = [15.143 0.998 1.4286 15.249 6.284 1.714 1.85 0.2953 27.657]; %L
Va = volumes(1)*0.3; %Arterial side
Vv = volumes(1)*0.7; %Venous side
Q1 = flowrates(1); V1 = volumes(2); %Lung
Q2 = flowrates(2); V2 = volumes(3); %Brain
Q3 = flowrates(3); V3 = volumes(4); %Visceral Fat
Q4 = flowrates(4); V4 = volumes(5); %SC Fat
Q5 = flowrates(5); V5 = volumes(6); %
Q6 = flowrates(6); V6 = volumes(7); %Lung
Q7 = flowrates(7); V7 = volumes(8); %Lung
Q8 = flowrates(8); V8 = volumes(9); %Lung

%The half life
t_half = parameters(1);
%The fraction unbound of drug
fu = parameters(4);
%Water-octanol partition
P = parameters(5);
F_brain = P/(1+P); %fraction into brain
%Amount eliminated by urine
F_urine = parameters(6);

tspan = [0 t_half*5];
y0 = [0 Dose/Vv 0 0 0 0 0 0 0 0]; 

x = [];
y = [];
z = [];
%x(1) = Arterial Blood Flow %x(2) = Venous Blood flow
%x(3) = Lung %x(4) = Brain
%x(5) = Visceral Fat %x(6) = SC Fat
%x(7) = Liver %x(8) = GI
%x(9) = Kidney %x(10) = Muscle
for ke_K = linspace(0.5, 10, 10)
    for ke_L = linspace(0.05, 10, 10)
        ode_MC = @(t,x) [(Q1/Va)*x(3) - x(1)*fu*(F_brain*Q2+Q3+Q4+Q5+Q6+Q7+Q8)/Va;... %Arterial
            -(Q1/Vv)*fu*x(2) + (1-F_brain)*(Q2/Vv)*x(4) + ((Q3+0.1*Q6)/Vv)*x(5) + (Q4/Vv)*x(6) + ((Q5+0.9*Q6)/Vv)*x(7) + (Q7/Vv)*x(9) + (Q8/Vv)*x(10);... %Venous
            (Q1/V1)*(fu*x(2) - x(3));... %Lung
            (Q2/V2)*(fu*F_brain*x(1) - (1-F_brain)*x(4));... %Brain
            (Q3/V3)*(fu*x(1) - x(5)) + 0.1*(Q6/V3)*(x(8)-x(5));... %Visceral
            (Q4/V4)*(fu*x(1) - x(6));... %SC Fat
            (Q5/V5)*(fu*x(1) - x(7)) + 0.9*(Q6/V5)*(x(8)-x(7)) - ke_L*fu*x(7);... %Liver
            (Q6/V6)*(fu*x(1) - x(8));... %GI
            (Q7/V7)*(fu*x(1) - x(9)) - ke_K*fu*x(9);... %Kindey
            (Q8/V8)*(fu*x(1) - x(10))]; %Muscle
        [t, C] = ode45(ode_MC, tspan, y0);
        logC = log(C);
        j = 1;
        while t(j) <  2*t_half
            j = j + 1;
        end %This searches for a timepoint deep into the profile to avoid distribution effects
        j2 = 1;
        while t(j2) < 3*t_half
            j2 = j2 + 1;
        end
        q = polyfit(t(j:j2),logC(j:j2, 2),1); %fits a line to the log plot slope
        slope = q(1);
        x = [x, ke_K];
        y = [y, ke_L];
        z = [z, slope];
    end
end

x = x';
y = y';
z = z';
coeff = coeffvalues(fit([x, y], z, 'poly22'));
disp(coeff)

est = coeff(1) + coeff(2)*x + coeff(3)*y + coeff(4)*x.^2 + coeff(6)*y.^2;
% figure(3)
% scatter3(x, y, z);
% hold on;
% scatter3(x, y, est);
% xlabel('ke_K');
% ylabel('ke_L');
% zlabel('slope');
% legend('Actual Slope', 'Predicted Slope', 'Location', 'best');
% hold off;

eq1 = @(k1) F_urine*(log(2)/t_half - coeff(1)) + coeff(2)*k1 + coeff(4)*k1^2;
eq2 = @(k2) (1-F_urine)*(log(2)/t_half - coeff(1)) + coeff(3)*k2 + coeff(6)*k2^2;

solvedKeK = fsolve(eq1, 2*log(2)/t_half);
solvedKeL = fsolve(eq2, 2*log(2)/t_half);
if solvedKeK < 0 || solvedKeL < 0
    error('Half life too big.')
end
end