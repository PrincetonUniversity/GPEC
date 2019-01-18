function w = getIntervalOptim(psi12, s1to3, t12)
    %These are the parameters I had been using for some time...
    gamma = .009;
    a_a = 30;
    a_s = 146;
    a_e = 3;
    b_a = 14316;
    b_s = 200000;
    b_e = 11000;

    w0 = [gamma; a_a; a_s; a_e; b_a; b_s; b_e];

    figure;
    subplot(2,1,1);
    plot(t12); hold on; plot(getT(w0));

    options = optimset('MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-9,'TolX',1e-9);
    w = fmincon(@esttime, w0,-eye(7),zeros(7,1),[],[],[],[],[],options);
    %w = fminsearch(@esttime, w0,options);  %fmincon, with constraints, was better!

    subplot(2,1,2);
    plot(t12); hold on; plot(getT(w));

    function rOut = esttime(y)
        tX = getT(y);
        rOut = norm(t12-tX);
        %rOut = norm(sign(t12).*t12.^2-sign(tX).*tX.^2); %Overemphasized peaks
    end

    function tOut = getT(wIn)
        gamma = wIn(1);
        a_a   = wIn(2);
        a_s   = wIn(3);
        a_e   = wIn(4);
        b_a   = wIn(5);
        b_s   = wIn(6);
        b_e   = wIn(7);

        tOut = gamma + a_a/b_a * abs(log(1+b_a*abs(psi12(:,1) - 0)) - ...
                                     log(1+b_a*abs(psi12(:,2) - 0))   );
        tOut = tOut + a_e/b_e * abs( ...
                        log(1+b_e*abs(psi12(:,1) - 1)) - ...
                        log(1+b_e*abs(psi12(:,2) - 1))   );
        for jS = 1:3
            tOut = tOut + a_s/b_s * abs( ...
                        log(1+b_s*abs(psi12(:,1) - s1to3(jS))) - ...
                        log(1+b_s*abs(psi12(:,2) - s1to3(jS)))   );
        end
    end
end