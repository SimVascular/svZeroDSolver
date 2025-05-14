clear all, close all, clc

syms r v S tau V Pin Pout Qin Qout c
syms r_ v_ S_ tau_ V_ Pin_ Pout_ Qin_ Qout_ c_
syms rho d Ro W1 W2 eta a a_plus sigma_o

  d1_state_vars=[Pin Qin Pout Qout r v S tau V];
  d2_state_vars=[Pin_ Qin_ Pout_ Qout_ r_ v_ S_ tau_ V_];
    C=@(r) ((1+(r/Ro))^2); 
    C_=@(r,r_) 2*(1+(r/Ro))*(1/Ro)*r_; 
            Res_pmp(1)=(rho*d*v_) + (d/Ro)*(1 + (r/Ro))*S - Pout*C(r); 
            Res_pmp(2)=-S + tau + 4*(1 - ((C(r))^-3))*(W1 + C(r)*W2) + 2*eta*C_(r,r_)*(1-2*((C(r))^-6));
            Res_pmp(3)=tau_ + (a*tau) - (sigma_o*a_plus);
            Res_pmp(4)=Qin - Qout - V_;
            
            Res_pmp(5)=4*pi*(Ro^2)*((1 + (r/Ro))^2)*v - V_;
            Res_pmp(6)=r_ - v;
            Res_pmp(7)=Pin - Pout;
            
            % display title for C equations C++ code
            fprintf('C++ code for C expressions\n\n');

            % display C equations C++ code
            for i = 1:length(Res_pmp)
                Res_pmp_text = ccode(Res_pmp(i));  % Convert the symbolic expression to a string of C++ executable code
                Res_pmp_text = strrep(Res_pmp_text, '3.141592653589793', 'M_PI');
                Res_pmp_text = regexprep(Res_pmp_text, '^\s*t\d+\s*=\s*', '');
                fprintf('system.C(global_eqn_ids[%d]) = %s\n', ...
                        i-1, Res_pmp_text);
            end

            % display title for dC_dy C++ code 
            fprintf(' \nC++ code for dC_dy expressions\n\n');

            % calculate jacobian
            for i=1:length(Res_pmp) %assembles jacobian matrix
                for j=1:length(d1_state_vars)
                    Jac1_pmp(i,j)=simplify(diff(Res_pmp(i),d1_state_vars(j)));
                    Jac1 = Jac1_pmp;

                    %display dC_dy C++ code
                    if Jac1(i,j) ~= 0
                        jacobian1_text = ccode(Jac1(i,j));  % Convert the symbolic expression to a string of C++ executable code
                        jacobian1_text = strrep(jacobian1_text, '3.141592653589793', 'M_PI');
                        jacobian1_text = regexprep(jacobian1_text, '^\s*t\d+\s*=\s*', '');
                        fprintf('system.dC_dy.coeffRef(global_eqn_ids[%d], global_var_ids[%d]) = %s\n', ...
                                i-1, j-1, jacobian1_text);
                    end
                end
            end
            
            % display title for dC_dydot C++ code
            fprintf(' \nC++ code for dC_dydot expressions\n\n');

            % calculate dC_dydot jacobian
            for i=1:length(Res_pmp) %assembles jacobian matrix
                for j=1:length(d2_state_vars)
                    Jac2_pmp(i,j) = simplify(diff(Res_pmp(i),d2_state_vars(j)));
                    Jac2 = Jac2_pmp;

                    % display d2C_dy2 C++ code
                    if Jac2(i,j) ~= 0
                        jacobian2_text = ccode(Jac2(i,j));  % Convert the symbolic expression to a string of C++ executable code
                        jacobian2_text = strrep(jacobian2_text, '3.141592653589793', 'M_PI');
                        jacobian2_text = regexprep(jacobian2_text, '^\s*t\d+\s*=\s*', '');
                        fprintf('system.dC_dydot.coeffRef(global_eqn_ids[%d], global_var_ids[%d]) = %s\n', ...
                                i-1, j-1, jacobian2_text);
                    end
                end
            end

    % un-comment to display original jacobian matrices 
    % Jac1 = Jac1_pmp
    % Jac2 = Jac2_pmp

    constants = [rho d Ro W1 W2 eta a sigma_o];

     % display title for defining constants C++ code
     fprintf('\nC++ code for defining constants\n\n');

    for i = 1:length(constants)
        fprintf('double %s = parameters[global_param_ids[ParamId::%s]];\n', constants(i), constants(i))
    end
    
     % display title for defining y variables C++ code
     fprintf('\nC++ code for defining y variables\n\n');

    for i = 1:length(d1_state_vars)
        fprintf('double %s = y[global_var_ids[%d]]; \n', d1_state_vars(i), i-1)
    end

    % display title for defining dy variables C++ code
     fprintf('\nC++ code for defining dy variables\n\n');

    for i = 1:length(d2_state_vars)
        fprintf('double %s = dy[global_var_ids[%d]]; \n', d2_state_vars(i), i-1)
    end
