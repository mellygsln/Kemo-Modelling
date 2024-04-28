function [t,L,S,K] = runge_kutta(t0,tN,Nt,l0,s0,k0,f,g,q)
    h = (tN-t0)/Nt;
    t = t0:h:tN;
    L(1) = l0;
    S(1) = s0;
    K(1) = k0;

    for j = 1:length(t)-1
        f1 = f(t(j), L(j), S(j), K(j));
        g1 = g(t(j), L(j), S(j), K(j));
        q1 = q(t(j), L(j), S(j), K(j));

        f2 = f(t(j)+0.5*h, L(j)+0.5*h*f1, S(j)+0.5*h*g1, K(j)+0.5*h*q1);
        g2 = g(t(j)+0.5*h, L(j)+0.5*h*f1, S(j)+0.5*h*g1, K(j)+0.5*h*q1);
        q2 = q(t(j)+0.5*h, L(j)+0.5*h*f1, S(j)+0.5*h*g1, K(j)+0.5*h*q1);

        f3 = f(t(j)+0.5*h, L(j)+0.5*h*f2, S(j)+0.5*h*g2, K(j)+0.5*h*q2);
        g3 = g(t(j)+0.5*h, L(j)+0.5*h*f2, S(j)+0.5*h*g2, K(j)+0.5*h*q2);
        q3 = q(t(j)+0.5*h, L(j)+0.5*h*f2, S(j)+0.5*h*g2, K(j)+0.5*h*q2);
        
        f4 = f(t(j)+h, L(j)+h*f3, S(j)+h*g3, K(j)+h*q3);
        g4 = g(t(j)+h, L(j)+h*f3, S(j)+h*g3, K(j)+h*q3);
        q4 = q(t(j)+h, L(j)+h*f3, S(j)+h*g3, K(j)+h*q3);
        L(j+1) = L(j) + (h/6)*(f1+2*f2+2*f3+f4);
        S(j+1) = S(j) + (h/6)*(g1+2*g2+2*g3+g4);
        K(j+1) = K(j) + (h/6)*(q1+2*q2+2*q3+q4);
    end
end