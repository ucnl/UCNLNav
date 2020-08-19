function [ x_best, y_best, rerr, it_cnt, f_ev_cnt ] = Nav_hjs_2d_solve_demo(eps, base_elements, x_prev, y_prev, z_prev,...
    max_iterations, precision_threshold, step_size)  

var_num = 2;
x0 = [ x_prev y_prev ];
f0 = eps(base_elements, x_prev, y_prev, z_prev);
deltas0 = [ step_size step_size ]; 
deltas = deltas0;
it_cnt = 0;
f_ev_cnt = 1;

while ((max(deltas) > precision_threshold) &&...
       (it_cnt < max_iterations))
   
       [x1, f1, deltas, ecnt] = hjs_explore(x0, f0, deltas, var_num);
       f_ev_cnt = f_ev_cnt + ecnt;
            
       while ((f1 < f0) && (it_cnt < max_iterations))
 
           deltas = deltas0;
           x2 = x1 * 2 - x0; 
           f2 = eps(base_elements, x2(1), x2(2), z_prev);
           f_ev_cnt = f_ev_cnt + 1;
           [x2, f2, dts, ecnt] = hjs_explore(x2, f2, deltas, var_num);
           f_ev_cnt = f_ev_cnt + ecnt;
           
           line([x0(1) x1(1) x2(1)], [x0(2) x1(2) x2(2)], 'color', 'green');
           pause(0.01);                      
           
           if (f2 > f1)                                           
              x0 = x1;
              f0 = f1;              
           else         
              x0 = x1;
              f0 = f1;
              x1 = x2;
              f1 = f2;
              
              plot(x1(1), x1(2), 'g+');
              it_cnt = it_cnt + 1;
           end
       end
   
       
   it_cnt = it_cnt + 1;
end

x_best = x0(1);
y_best = x0(2);
rerr = sqrt(f0);

    function [x_p, f_p, dlts2, fecnt] = hjs_explore(x_c, f_c, dlts, v_num)
       
       x_p = x_c;
       f_p = f_c;
       dlts2 = dlts;
       fecnt = 0;
       
       plot(x_c(1), x_c(2), '.g');
       
       for nn = 1:v_num
           
           x_c_prev = x_p(nn);
           f_prev = f_p;
           
           x_p(nn) = x_c_prev + dlts(nn);
           f_p = eps(base_elements, x_p(1), x_p(2), z_prev);
           fecnt = fecnt + 1;
                      
           if (f_p > f_c)
               
               f_p = f_prev;
               x_p(nn) = x_c_prev;
               
               plot(x_p(1), x_p(2), '.g');
                              
               x_p(nn) = x_c_prev - dlts(nn);
               f_p = eps(base_elements, x_p(1), x_p(2), z_prev);   
               fecnt = fecnt + 1;

               if (f_p > f_c)
                    
                   f_p = f_prev;
                   x_p(nn) = x_c_prev;
                   
                   plot(x_p(1), x_p(2), '.g');
                   dlts2(nn) = dlts2(nn) / 2;
                   
               else
                                      
                   plot(x_p(1), x_p(2), '+g');
                   
               end     
               
           else
               
               plot(x_p(1), x_p(2), '+g');             
               
           end     
       end               
    end

end