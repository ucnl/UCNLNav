function [ x_best, y_best, z_best, rerr, it_cnt ] = Nav_nlm_3d_solve(eps, base_elements, x_prev, y_prev, z_prev,...
    max_iterations, precision_threshold, simplex_size)  
   
    NLM_A = 1.0;
    NLM_B = 0.5;
    NLM_R = 0.5;
    NLM_Q = 0.5;
    NLM_G = 2.0;

    xix = zeros(1, 4);
    xiy = zeros(1, 4);
    xia = zeros(1, 4);
    fxi = zeros(1, 4);

    xix(1) = x_prev;
    xiy(1) = y_prev;
    xia(1) = z_prev;
    xix(2) = xix(1) + simplex_size;
    xiy(2) = xiy(1) + simplex_size;
    xia(2) = xia(1) - simplex_size;
    xix(3) = xix(1) - simplex_size / 2.0;
    xiy(3) = xiy(1) + simplex_size / 2.0;
    xia(3) = xia(1) + simplex_size / 2.0;
    xix(4) = xix(1) + simplex_size / 2.0;
    xiy(4) = xiy(1) - simplex_size / 2.0;
    xia(4) = xia(1) - simplex_size / 2.0;

    is_finished = false;
    it_cnt = 0;

    while ~is_finished
        
        for n = 1:4
            fxi(n) = eps(base_elements, xix(n), xiy(n), xia(n));
        end
        
        sigma = std(fxi);
        
        it_cnt = it_cnt + 1;
        if (it_cnt > max_iterations) || (sigma < precision_threshold)
           is_finished = true;
        else        
            % Sort vertices
            if fxi(1) > fxi(2) 
                tmp = fxi(1); fxi(1) = fxi(2); fxi(2) = tmp;
                tmp = xix(1); xix(1) = xix(2); xix(2) = tmp;
                tmp = xiy(1); xiy(1) = xiy(2); xiy(2) = tmp;
                tmp = xia(1); xia(1) = xia(2); xia(2) = tmp;
            end
            
            if fxi(1) > fxi(3)
                tmp = fxi(1); fxi(1) = fxi(3); fxi(3) = tmp;
                tmp = xix(1); xix(1) = xix(3); xix(3) = tmp;
                tmp = xiy(1); xiy(1) = xiy(3); xiy(3) = tmp;
                tmp = xia(1); xia(1) = xia(3); xia(3) = tmp;
            end
            
            if fxi(1) > fxi(4)
                tmp = fxi(1); fxi(1) = fxi(4); fxi(4) = tmp;
                tmp = xix(1); xix(1) = xix(4); xix(4) = tmp;
                tmp = xiy(1); xiy(1) = xiy(4); xiy(4) = tmp;
                tmp = xia(1); xia(1) = xia(4); xia(4) = tmp;
            end
            
            if fxi(2) > fxi(3)
                tmp = fxi(2); fxi(2) = fxi(3); fxi(3) = tmp;
                tmp = xix(2); xix(2) = xix(3); xix(3) = tmp;
                tmp = xiy(2); xiy(2) = xiy(3); xiy(3) = tmp;
                tmp = xia(2); xia(2) = xia(3); xia(3) = tmp;
            end
            
            if fxi(2) > fxi(4)
                tmp = fxi(2); fxi(2) = fxi(4); fxi(4) = tmp;
                tmp = xix(2); xix(2) = xix(4); xix(4) = tmp;
                tmp = xiy(2); xiy(2) = xiy(4); xiy(4) = tmp;
                tmp = xia(2); xia(2) = xia(4); xia(4) = tmp;
            end
            
            if fxi(3) > fxi(4)
                tmp = fxi(3); fxi(3) = fxi(4); fxi(4) = tmp;
                tmp = xix(3); xix(3) = xix(4); xix(4) = tmp;
                tmp = xiy(3); xiy(3) = xiy(4); xiy(4) = tmp;
                tmp = xia(3); xia(3) = xia(4); xia(4) = tmp;
            end

            % (2) Calculate x0, the centroid of all points except xn+1
            x0x = mean(xix(1:3));
            x0y = mean(xiy(1:3));
            x0a = mean(xia(1:3));

            % (3) reflect the xh (xi(4)) point to xr
            xrx = x0x + NLM_A * (x0x - xix(4));
            xry = x0y + NLM_A * (x0y - xiy(4));
            xra = x0a + NLM_A * (x0a - xia(4));

            % function value in the xr point                    
            fr = eps(base_elements, xrx, xry, xra);

            % if fx1 <= fxr <= fxn replace worst point xn+1 with xr
            if (fr >= fxi(1)) && (fxi(3) >= fr)
                xix(4) = xrx;
                xiy(4) = xry;
                xia(4) = xra;
                fxi(4) = fr;
            else
                % (4) expansion
                if fr < fxi(1)
                    xex = x0x + NLM_G * (xrx - x0x);
                    xey = x0x + NLM_G * (xry - x0x);
                    xea = x0x + NLM_G * (xra - x0x);
                    fe = eps(base_elements, xex, xey, xea);

                    if fe < fr
                        xix(4) = xex;
                        xiy(4) = xey;
                        xia(4) = xea;
                        fxi(4) = fe;
                        % to step 1
                    else
                        xix(4) = xrx;
                        xiy(4) = xry;
                        xia(4) = xra;
                        fxi(4) = fr;
                        % to step 1
                    end                
                else
                    xcx = x0x + NLM_R * (xix(4) - x0x);
                    xcy = x0y + NLM_R * (xiy(4) - x0y);
                    xca = x0a + NLM_R * (xia(4) - x0a);
                    fc = eps(base_elements, xcx, xcy, xca);

                    if fc < fxi(4)
                        xix(4) = xcx;
                        xiy(4) = xcy;
                        xia(4) = xca;
                        fxi(4) = fc;
                        % to step 1
                    else
                        xix(2) = xix(1) + NLM_Q * (xix(2) - xix(1));
                        xiy(2) = xiy(1) + NLM_Q * (xiy(2) - xiy(1));
                        xia(2) = xia(1) + NLM_Q * (xia(2) - xia(1));

                        xix(3) = xix(1) + NLM_Q * (xix(3) - xix(1));
                        xiy(3) = xiy(1) + NLM_Q * (xiy(3) - xiy(1));
                        xia(3) = xia(1) + NLM_Q * (xia(3) - xia(1));

                        xix(4) = xix(1) + NLM_Q * (xix(4) - xix(1));
                        xiy(4) = xiy(1) + NLM_Q * (xiy(4) - xiy(1));
                        xia(4) = xia(1) + NLM_Q * (xia(4) - xia(1));
                    end
                end
            end
        end
    end
    
    x_best = xix(1);
    y_best = xiy(1);
    z_best = xia(1);
    rerr = sqrt(eps(base_elements, xix(1), xiy(1), xia(1)));
    
end