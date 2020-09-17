
function [] = enzyme_kinetics(t, x)
% Enzyme kinetics - the Michaelis Menten Model for reversible enzyme-substrate and enzyme-product reactions
    function [f] = dXdT(t, x)
        S = x(1);
        P = x(2);
        E = x(3);
        ES = x(4);

        k1 = 10;
        k2 = 1;
        km2 = 0.5;
        km1 = 1;

        dSdt = -k1*E*S + km1*ES;
        dPdt = k2*ES - km2*E*P;
        dEdt = -k1*E*S + km1*ES + k2*ES - km2*E*P;
        dESdt = k1*E*S +km2*E*P - km1*ES - k2*ES;
        f = [dSdt; dPdt; dEdt; dESdt]; 
    end
[t,y] = ode15s(@dXdT,t, x);
figure
plot(t, y,'linewidth', 2)
legend({'[S]', '[P]', '[E]', '[ES]'},'Location','bestoutside')
legend('boxoff')
title("Enzyme Kinetics - the Michaelis Menten Model")
xlabel("t")
ylabel("[C]")
end