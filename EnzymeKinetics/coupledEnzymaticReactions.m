
function [] = coupledEnzymaticReactions(t, x)
% Only step 1 and step 6 of glycolysis is studied here.
% Take the following initial values -  [Glu]0 = 12.87, [Glu-6-P = Y]0 = 1.5, [ATP]= 2.5, [ADP] = 1.2
    function [f] = dXdT(t, x)
        X = x(1);
        Y = x(2);
        ATP = x(3);
        ADP = x(4);

        vmax1 = 1398;
        kATP1 = 0.1;
        kGlu1 = 0.37;
        k6 = 68.48;
        
        rs1 = vmax1*X*ATP/(1 + ATP/kATP1 + X/kGlu1);
        rs6 = k6*ADP;

        dXdt = -rs1;
        dYdt = rs1;
        dATPdt = -rs1 + rs6;
        dADPdt = rs1 - rs6;
        f = [dXdt; dYdt; dATPdt; dADPdt]; 
    end
[t,y] = ode15s(@dXdT,t, x);
figure
plot(t, y,'linewidth', 2)
legend({'[Glu]', '[Glu-6-P]', '[ATP]', '[ADP]'},'Location','bestoutside')
legend('boxoff')
title("Coupled Enzymatic Reactions")
xlabel("t (min)")
ylabel("Concentrations")
end
