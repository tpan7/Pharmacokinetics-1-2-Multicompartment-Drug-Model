function kA = ka_solveMC(parameters, Dose, ke_K, ke_L, route)
%1 for Intravenous, 2 for Oral, 3 for Intramusclar, 4 for Subcutaneous
%5 for Buccal, 6 for Sublingual, 7 for Lung, 8 for Rectal, 9 for
%transdermal

%The half life
t_half = parameters(1);
%The bioavailability factor
F = parameters(2);
%Time where max concentration occurs
t_max = parameters(3);
%The fraction unbound of drug
fu = parameters(4);
%Water-octanol partition
P = parameters(5);
F_brain = P/(1+P); %fraction into brain

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

tol = .1;
low_guess_ka = 0; %define some high and low bounds for guessing
high_guess_ka = 50;

T_max = t_max/2; %define some arbitrary T_max for loop to function
tspan = [0 t_half*5];
y0 = [0 0 0 0 0 0 0 0 0 0 Dose]; 
while T_max < t_max - tol * t_max || T_max > t_max + tol * t_max
    ka = (low_guess_ka + high_guess_ka)/2;
    if ka == 0 || ka == 50
        error('Converged to the high or low bound. Needs to be changed in ka_solveMC.m')
    end
    ode_MC = @(t,x) [(Q1/Va)*x(3) - x(1)*fu*(F_brain*Q2+Q3+Q4+Q5+Q6+Q7+Q8)/Va;... %Arterial
            -(Q1/Vv)*fu*x(2)+(1-F_brain)*(Q2/Vv)*x(4)+((Q3+0.1*Q6)/Vv)*x(5)+(Q4/Vv)*x(6)+((Q5+0.9*Q6)/Vv)*x(7)+(Q7/Vv)*x(9)+(Q8/Vv)*x(10);... %Venous
            F*ka*x(11)*route(7)/V1 + (Q1/V1)*(fu*x(2) - x(3));... %Lung
            (Q2/V2)*(fu*F_brain*x(1) - (1-F_brain)*x(4));... %Brain
            (Q3/V3)*(fu*x(1) - x(5)) + 0.1*(Q6/V3)*(x(8)-x(5));... %Visceral
            F*ka*x(11)*route(9)/V4 + (Q4/V4)*(fu*x(1) - x(6));... %SC Fat
            (Q5/V5)*(fu*x(1) - x(7)) + 0.9*(Q6/V5)*(x(8)-x(7)) - ke_L*fu*x(7);... %Liver
            F*ka*x(11)*(route(2)+route(5)+route(6)+route(8))/V6 + (Q6/V6)*(fu*x(1) - x(8));... %GI
            (Q7/V7)*(fu*x(1) - x(9)) - ke_K*fu*x(9);... %Kindey
            (Q8/V8)*(fu*x(1) - x(10));... %Muscle
            -ka*x(11)]; %Balance for the GI tract
    [t, C] = ode45(ode_MC, tspan, y0);
    [M, I] = max(C(:,2));
    T_max = t(I);
    if T_max < t_max - tol * t_max
        high_guess_ka = ka;
    elseif T_max > t_max + tol * t_max
        low_guess_ka = ka;
    end
end
kA = ka;
end