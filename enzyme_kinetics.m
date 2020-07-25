
function [] = enzyme_kinetics(t, x)
% Enzyme kinetics - the Michaelis Menteen Model
    function [f] = dXdT(t, x)
        S = x(1);
        P = x(2);
        E = x(3);
        ES = x(4);

        k1 = 20;
        k2 = 1;
        km1 = 1;

        dSdt = -k1*E*S + km1*ES;
        dPdt = k2*ES;
        dEdt = -k1*E*S + km1*ES + k2*ES;
        dESdt = k1*E*S - km1*ES - k2*ES;
        f = [dSdt; dPdt; dEdt; dESdt]; 
    end
[t,y] = ode15s(@dXdT,t, x);
figure
plot(t, y,'linewidth', 2)
legend({'[S]', '[P]', '[E]', '[ES]'},'Location','bestoutside')
legend('boxoff')
title("Enzyme Kinetics - the Michaelis Menteen Model")
xlabel("t")
ylabel("[C]")
end