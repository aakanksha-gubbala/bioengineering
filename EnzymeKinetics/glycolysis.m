
function [] = glycolysis(t, x)
% X = Glu and Y = Glu-6-P and Z = Fructo-6-P and W = Fructo-1,6-P
% Take the following initial concentrations:
% [Glu]0 = 12.874, [Glu-6-P]0 = 1
% [Fructo-6-P]0 = 0, [Fructo-1, 6-P]0 = 0
% [P1]0 = 0, [P2]0 = 0
% [ATP]0 = 2.1, [ADP]0 = 1.4, [AMP]0 = 0
    function [f] = dGlydT(t, x)
        X = x(1);
        Y = x(2);
        Z = x(3);
        W = x(4);
        ATP = x(5);
        ADP = x(6);
        AMP = x(7);
        P1 = x(8);
        P2 = x(9);

        vmax1 = 1398;
        kATP1 = 0.1;
        kGlu1 = 0.37;
        
        kGlu6P3 = 0.80;
        kFru6P3 = 0.15;
        vfmax3 = 140.482;
        vbmax3 = 140.282;
        
        kFruc6P4 = 1;
        K = 1;
        vmax4 = 1;
        
        k2 = 2.26;
        k5 = 6.04662;
        k6 = 68.48;
        k7 = 3.21;
        k8f = 432.8;
        k8b = 1;
        
        
        rs1 = vmax1*X*ATP/(1 + ATP/kATP1 + X/kGlu1);
        rs2 = k2*ATP*Y;
        rs3 = (vfmax3*Y/kGlu6P3 - vbmax3*Z/kFru6P3)/(1 + Y/kGlu6P3 + Z/kFru6P3);
        rs4 = vmax4*Y^2/(kFruc6P4*(1 + K*(ATP/AMP)^2 + Z^2));
        rs5 = k5*W;
        rs6 = k6*ADP;
        rs7 = k7*ATP;
        rs8 = k8f*ATP*AMP - k8b*ADP^2;

        dXdt = -rs1;
        dYdt = rs1 - rs2 - rs3;
        dZdt = rs3 -rs4;
        dWdt = rs4 - rs5;
        dATPdt = -rs1 - rs2 - rs4 + rs6 - rs7 - rs8;
        dADPdt = rs1 + rs2 + rs4 - rs6 + rs7 + 2*rs8;
        dAMPdt = -rs8;
        dP1dt = rs2;
        dP2dt = rs5;
        f = [dXdt; dYdt; dZdt; dWdt; dATPdt; dADPdt; dAMPdt; dP1dt; dP2dt]; 
    end
[t,y] = ode15s(@dGlydT,t, x);
figure
plot(t, y,'linewidth', 1)
legend({'[Glu]', '[Glu-6-P]', '[Fructo-6-P]','[Fructo-1,6-P]','[ATP]', '[ADP]', '[AMP]','[$P_1$]', '[$P_2$]'},'interpreter', 'latex','Location','bestoutside')
legend('boxoff')
title("Glycolysis Reactions")
xlabel("t (min)")
ylabel("Concentration (mM)")
end
