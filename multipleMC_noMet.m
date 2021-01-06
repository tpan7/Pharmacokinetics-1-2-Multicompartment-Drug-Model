function multipleMC_noMet(parameters, route, num_doses)
%Which route of administration
%1 for Intravenous, 2 for Oral, 3 for Intramusclar, 4 for Subcutaneous
%5 for Buccal, 6 for Sublingual, 7 for Lung, 8 for Rectal, 9 for
%transdermal
[temp, indx] = max(route);

%x(1) = Arterial Blood Flow %x(2) = Venous Blood flow
%x(3) = Lung %x(4) = Brain
%x(5) = Visceral Fat %x(6) = SC Fat
%x(7) = Liver %x(8) = GI
%x(9) = Kidney %x(10) = Muscle
flowrates = 60*[4.7 0.750 0.3413 0.0597 0.3917 0.7838 1 0.705]; %L/hr
volumes = [15.143 0.998 1.4286 15.249 6.284 1.714 1.85 0.2953 27.657]; %L

Va = volumes(1)*0.3; %Arterial side
Vv = volumes(1)*0.7; %Venous side
Q1 = flowrates(1); V1 = volumes(2);
Q2 = flowrates(2); V2 = volumes(3);
Q3 = flowrates(3); V3 = volumes(4);
Q4 = flowrates(4); V4 = volumes(5);
Q5 = flowrates(5); V5 = volumes(6);
Q6 = flowrates(6); V6 = volumes(7);
Q7 = flowrates(7); V7 = volumes(8);
Q8 = flowrates(8); V8 = volumes(9);

%The bioavailability factor
F = parameters(2);
if indx == 1
    F = 1;
end
%The fraction unbound of drug
fu = parameters(4);
%Water-octanol partition
P = parameters(5);
%Fraction that goes into brain
F_brain = P/(1+P);
%The end of the plotting
range = parameters(7);

dosing_times = [0];
Dose = [];
for i=1:num_doses-1
    ti = input(['At what time did you give dose #' num2str(i+1) ' (hr)? (Assuming dose #1 was at t=0): ']);
    if ti > range
        error('This dose is outside the time you are integrating over.')
    elseif ti <= dosing_times(end)
        error('This dose time is before or at the same time as the previous one.'); 
    end
    dosing_times = [dosing_times ti];
end
for i=1:num_doses
    Dosei = input(['What was dose #' num2str(i) ' (mg)?: ']);
    Dose = [Dose Dosei];
end

%Get the elimination terms
[ke_K, ke_L] = ke_solveMC(parameters, Dose(1));

if indx ~= 1 && indx ~= 3 && indx ~= 4
    ka = ka_solveMC(parameters, Dose(1), ke_K, ke_L, route);
end

if indx == 1 || indx == 3 || indx == 4
    t = [];
    C = [[]];
    for i=1:num_doses-1
        sys = @(t, x) [(Q1/Va)*x(3) - x(1)*fu*(F_brain*Q2+Q3+Q4+Q5+Q6+Q7+Q8)/Va;... %Arterial
            -(Q1/Vv)*fu*x(2) + (1-F_brain)*(Q2/Vv)*x(4) + ((Q3+0.1*Q6)/Vv)*x(5) + (Q4/Vv)*x(6) + ((Q5+0.9*Q6)/Vv)*x(7) + (Q7/Vv)*x(9) + (Q8/Vv)*x(10);... %Venous
            (Q1/V1)*(fu*x(2) - x(3));... %Lung
            (Q2/V2)*(fu*F_brain*x(1) - (1-F_brain)*x(4));... %Brain
            (Q3/V3)*(fu*x(1) - x(5)) + 0.1*(Q6/V3)*(x(8)-x(5));... %Visceral
            (Q4/V4)*(fu*x(1) - x(6));... %SC Fat
            (Q5/V5)*(fu*x(1) - x(7)) + 0.9*(Q6/V5)*(x(8)-x(7)) - ke_L*fu*x(7);... %Liver
            (Q6/V6)*(fu*x(1) - x(8));... %GI
            (Q7/V7)*(fu*x(1) - x(9)) - ke_K*fu*x(9);... %Kindey
            (Q8/V8)*(fu*x(1) - x(10))]; %Muscle
        tspan = [dosing_times(i) dosing_times(i+1)];
        if i == 1
            y = [0 route(1)*Dose(i)/Vv 0 0 0 route(4)*Dose(i)/V4 0 0 0 route(3)*Dose(i)/V8];
        else
            y = [C(end,1) C(end,2)+route(1)*Dose(i)/Vv C(end,3) C(end,4) ...
                C(end,5) C(end,6)+route(4)*Dose(i)/V4 C(end,7) C(end,8) C(end,9) C(end,10)+route(3)*Dose(i)/V8];
        end
        [temp_t, temp_C] = ode45(sys, tspan, y);
        t = [t; temp_t];
        C = [C; temp_C];
    end
    sys = @(t, x) [(Q1/Va)*x(3) - x(1)*fu*(F_brain*Q2+Q3+Q4+Q5+Q6+Q7+Q8)/Va;... %Arterial
            -(Q1/Vv)*fu*x(2) + (1-F_brain)*(Q2/Vv)*x(4) + ((Q3+0.1*Q6)/Vv)*x(5) + (Q4/Vv)*x(6) + ((Q5+0.9*Q6)/Vv)*x(7) + (Q7/Vv)*x(9) + (Q8/Vv)*x(10);... %Venous
            (Q1/V1)*(fu*x(2) - x(3));... %Lung
            (Q2/V2)*(fu*F_brain*x(1) - (1-F_brain)*x(4));... %Brain
            (Q3/V3)*(fu*x(1) - x(5)) + 0.1*(Q6/V3)*(x(8)-x(5));... %Visceral
            (Q4/V4)*(fu*x(1) - x(6));... %SC Fat
            (Q5/V5)*(fu*x(1) - x(7)) + 0.9*(Q6/V5)*(x(8)-x(7)) - ke_L*fu*x(7);... %Liver
            (Q6/V6)*(fu*x(1) - x(8));... %GI
            (Q7/V7)*(fu*x(1) - x(9)) - ke_K*fu*x(9);... %Kindey
            (Q8/V8)*(fu*x(1) - x(10))]; %Muscle
    tspan = [dosing_times(num_doses) range];
    y = [C(end,1) C(end,2)+route(1)*Dose(num_doses)/Vv C(end,3) C(end,4) ...
         C(end,5) C(end,6)+route(4)*Dose(num_doses)/V4 C(end,7) C(end,8) C(end,9) ...
         C(end,10)+route(3)*Dose(num_doses)/V8];
    [temp_t, temp_C] = ode45(sys, tspan, y);
    t = [t; temp_t];
    C = [C; temp_C];
else 
    t = [];
    C = [[]];
    for i=1:num_doses-1
        sys = @(t, x) [(Q1/Va)*x(3) - x(1)*fu*(F_brain*Q2+Q3+Q4+Q5+Q6+Q7+Q8)/Va;... %Arterial
            -(Q1/Vv)*fu*x(2)+(1-F_brain)*(Q2/Vv)*x(4)+((Q3+0.1*Q6)/Vv)*x(5)+(Q4/Vv)*x(6)+((Q5+0.9*Q6)/Vv)*x(7)+(Q7/Vv)*x(9)+(Q8/Vv)*x(10);... %Venous
            F*ka*x(11)*route(7)/V1 + (Q1/V1)*(fu*x(2) - x(3));... %Lung
            (Q2/V2)*(fu*F_brain*x(1) - (1-F_brain)*x(4));... %Brain
            (Q3/V3)*(fu*x(1) - x(5)) + 0.1*(Q6/V3)*(x(8)-x(5));... %Visceral
            F*ka*x(11)*route(9)/V4 + (Q4/V4)*(fu*x(1) - x(6));... %SC Fat
            (Q5/V5)*(fu*x(1) - x(7)) + 0.9*(Q6/V5)*(x(8)-x(7)) - ke_L*fu*x(7);... %Liver
            F*ka*x(11)*(route(2)+route(5)+route(6)+route(8))/V6 + (Q6/V6)*(fu*x(1) - x(8));... %GI
            (Q7/V7)*(fu*x(1) - x(9)) - ke_K*fu*x(9);... %Kindey
            (Q8/V8)*(fu*x(1) - x(10));... %Muscle
            -ka*x(11)]; %Balance for external
        tspan = [dosing_times(i) dosing_times(i+1)];
        if i == 1
            y = [0 0 0 0 0 0 0 0 0 0 Dose(i)];
        else
            y = [C(end,1) C(end,2) C(end,3) C(end,4) C(end,5) C(end,6) C(end,7) C(end,8) C(end,9) C(end,10) C(end,11)+Dose(i)];
        end
        [temp_t, temp_C] = ode45(sys, tspan, y);
        t = [t; temp_t];
        C = [C; temp_C];
    end
    sys = @(t, x) [(Q1/Va)*x(3) - x(1)*fu*(F_brain*Q2+Q3+Q4+Q5+Q6+Q7+Q8)/Va;... %Arterial
            -(Q1/Vv)*fu*x(2)+(1-F_brain)*(Q2/Vv)*x(4)+((Q3+0.1*Q6)/Vv)*x(5)+(Q4/Vv)*x(6)+((Q5+0.9*Q6)/Vv)*x(7)+(Q7/Vv)*x(9)+(Q8/Vv)*x(10);... %Venous
            F*ka*x(11)*route(7)/V1 + (Q1/V1)*(fu*x(2) - x(3));... %Lung
            (Q2/V2)*(fu*F_brain*x(1) - (1-F_brain)*x(4));... %Brain
            (Q3/V3)*(fu*x(1) - x(5)) + 0.1*(Q6/V3)*(x(8)-x(5));... %Visceral
            F*ka*x(11)*route(9)/V4 + (Q4/V4)*(fu*x(1) - x(6));... %SC Fat
            (Q5/V5)*(fu*x(1) - x(7)) + 0.9*(Q6/V5)*(x(8)-x(7)) - ke_L*fu*x(7);... %Liver
            F*ka*x(11)*(route(2)+route(5)+route(6)+route(8))/V6 + (Q6/V6)*(fu*x(1) - x(8));... %GI
            (Q7/V7)*(fu*x(1) - x(9)) - ke_K*fu*x(9);... %Kindey
            (Q8/V8)*(fu*x(1) - x(10));... %Muscle
            -ka*x(11)]; %Balance for external
    tspan = [dosing_times(num_doses) range];
    y = [C(end,1) C(end,2) C(end,3) C(end,4) C(end,5) C(end,6) C(end,7) C(end,8) C(end,9) C(end,10) C(end,11)+Dose(num_doses)];
    [temp_t, temp_C] = ode45(sys, tspan, y);
    t = [t; temp_t];
    C = [C; temp_C];
end

figure(1)
plot(t, C(:,1), 'b');
hold on;
plot(t, C(:,2), 'r');
plot(t, C(:,3), 'g');
plot(t, C(:,4), 'k');
plot(t, C(:,5), 'm');
plot(t, C(:,6), 'b--');
plot(t, C(:,7), 'r--');
plot(t, C(:,8), 'g--');
plot(t, C(:,9), 'k--');
plot(t, C(:,10), 'm--');
grid on;
xlabel('Time (hr)')
ylabel('Conc of drug(mg/L)');
switch indx
    case 1
        title('Drug Conc vs Time from IV Administration')
    case 2
        title('Drug Conc vs Time from Oral Administration')
    case 3
        title('Drug Conc vs Time from Intramusclar Administration')
    case 4
        title('Drug Conc vs Time from Subcutaneous Administration')
    case 5
        title('Drug Conc vs Time from Buccal Administration')
    case 6
        title('Drug Conc vs Time from Sublingual Administration')
    case 7
        title('Drug Conc vs Time from Lung Administration')
    case 8
        title('Drug Conc vs Time from Rectal Administration')
end
legend('Arterial Flow', 'Venous Flow', 'Lung', 'Brain', 'Visceral Fat', 'Subcutaneous Fat', 'Liver', 'GI Organs', 'Kidney', 'Muscle', 'Location', 'best')
xlim([0 range])
hold off

figure(2)
plot(t, log(C(:,1)), 'b')
hold on; grid on;
plot(t, log(C(:,2)), 'r')
plot(t, log(C(:,3)), 'g')
plot(t, log(C(:,4)), 'k')
plot(t, log(C(:,5)), 'm')
plot(t, log(C(:,6)), 'b--')
plot(t, log(C(:,7)), 'r--')
plot(t, log(C(:,8)), 'g--')
plot(t, log(C(:,9)), 'k--')
plot(t, log(C(:,10)), 'm--')
xlabel('Time (hr)')
ylabel('Log conc (mg/L)');
switch indx
    case 1
        title('Drug Log Conc vs Time from IV Administration')
    case 2
        title('Drug Log Conc vs Time from Oral Administration')
    case 3
        title('Drug Log Conc vs Time from Intramusclar Administration')
    case 4
        title('Drug Log Conc vs Time from Subcutaneous Administration')
    case 5
        title('Drug Log Conc vs Time from Buccal Administration')
    case 6
        title('Drug Log Conc vs Time from Sublingual Administration')
    case 7
        title('Drug Log Conc vs Time from Lung Administration')
    case 8
        title('Drug Log Conc vs Time from Rectal Administration')
end
legend('Arterial Flow', 'Venous Flow', 'Lung', 'Brain', 'Visceral Fat', 'Subcutaneous Fat', 'Liver', 'GI Organs', 'Kidney', 'Muscle', 'Location', 'best')
xlim([0 range])
ylim([-7 3])
hold off
end