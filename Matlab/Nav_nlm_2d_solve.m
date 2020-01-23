function [ x_best, y_best, rerr, it_cnt ] = Nav_nlm_2d_solve(eps, base_elements, x_prev, y_prev, z,...
                  max_iterations, precision_threshold, simplex_size)

NLM_A = 1.0;
NLM_B = 0.5;
NLM_G = 2.0;              
              
is_finished = false;

xix = zeros(1, 3);
xiy = zeros(1, 3);
fxi = zeros(1, 3);    

it_cnt = 0;

xix(1) = x_prev;
xiy(1) = y_prev;
xix(2) = xix(1) + simplex_size;
xiy(2) = xiy(1) + simplex_size;
xix(3) = xix(1) - simplex_size / 2.0;
xiy(3) = xiy(1) + simplex_size / 2.0;

while ~is_finished 

    fxi(1) = eps(base_elements, xix(1), xiy(1), z);
    fxi(2) = eps(base_elements, xix(2), xiy(2), z);
    fxi(3) = eps(base_elements, xix(3), xiy(3), z);

    if fxi(1) > fxi(2)
        tmp = fxi(1); fxi(1) = fxi(2); fxi(2) = tmp;
        tmp = xix(1); xix(1) = xix(2); xix(2) = tmp;
        tmp = xiy(1); xiy(1) = xiy(2); xiy(2) = tmp;
    end

    if fxi(1) > fxi(3)
        tmp = fxi(1); fxi(1) = fxi(3); fxi(3) = tmp;
        tmp = xix(1); xix(1) = xix(3); xix(3) = tmp;
        tmp = xiy(1); xiy(1) = xiy(3); xiy(3) = tmp;
    end

    if fxi(2) > fxi(3)
        tmp = fxi(2); fxi(2) = fxi(3); fxi(3) = tmp;
        tmp = xix(2); xix(2) = xix(3); xix(3) = tmp;
        tmp = xiy(2); xiy(2) = xiy(3); xiy(3) = tmp;
    end

    fl = fxi(1);
    fg = fxi(2);
    fh = fxi(3);

    xcx = (xix(1) + xix(2)) / 2.0;
    xcy = (xiy(1) + xiy(2)) / 2.0;

    xrx = (1.0 + NLM_A) * xcx - NLM_A * xix(3);
    xry = (1.0 + NLM_A) * xcy - NLM_A * xiy(3);

    fr = eps(base_elements, xrx, xry, z);

    if fr < fl 
        xex = (1.0 - NLM_G) * xcx + NLM_G * xrx;
        xey = (1.0 - NLM_G) * xcy + NLM_G * xry;

        fe = eps(base_elements, xex, xey, z);

        if fe < fr 
            xix(3) = xex;
            xiy(3) = xey;
        else
            xix(3) = xrx;
            xiy(3) = xry;
        end
    else        
        if (fr > fl) && (fr < fg) 
            xix(3) = xrx;
            xiy(3) = xry;
        else
            if (fr > fg) && (fr < fh)
                xix(3) = xrx; % tmp = xix(3); xix(3) = xrx; xrx = tmp;
                xiy(3) = xry; % tmp = xiy(3); xiy(3) = xry; xry = tmp;
                fxi(3) = fr;  % tmp = fxi(3); fxi(3) = fr; fr = tmp;
            else
                if fh < fr
                    %
                end
            end

            xsx = NLM_B * xix(3) + (1.0 - NLM_B) * xcx;
            xsy = NLM_B * xiy(3) + (1.0 - NLM_B) * xcy;
            fs = eps(base_elements, xsx, xsy, z);

            if fs < fh
                xix(3) = xsx;
                xiy(3) = xsy;                
            else
                xix(2) = (xix(2) - xix(1)) / 2.0;
                xiy(2) = (xiy(2) - xiy(1)) / 2.0;
                xix(3) = (xix(3) - xix(1)) / 2.0;
                xiy(3) = (xiy(3) - xiy(1)) / 2.0;
            end
        end
    end

    it_cnt = it_cnt + 1;
    is_finished = (it_cnt >= max_iterations) || (std(fxi) <= precision_threshold);
end

x_best = xix(1);
y_best = xiy(1);
rerr = sqrt(eps(base_elements, xix(1), xiy(1), z));

end
