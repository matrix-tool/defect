function d = defect_h(U, METHOD, SV_TOLERANCE)
% 20170209 WB
% 20180824 W. Bruzda, name[at]uj.edu.pl : name = w.bruzda
%
% http://chaos.if.uj.edu.pl/~karol/hadamard/
% https://github.com/matrix-toolbox/
%
% [Restricted] defect of a unitary matrix U subjected to additional constraints.
%
% >> version % 9.1.0.441655 (R2016b)

    N = size(U, 1);
    if VERBOSE_MODE
        disp('Wait... Preparing matrix ''R''...');
    end
    z = 0; % number of zeros above the diagonal in U
    for j = 1 : N - 1
        for k = j + 1 : N
            z = z + (abs(U(j, k)) < 1e-12);
        end
    end
    % z = MUBs * N * (N - 1) / 2 for MUBs

    tau = N * (N - 1) / 2;
    trivial_phases = N - 1;     % equivalence relation...
    R = zeros(tau * 2, tau);    % system of linear equations matrix split onto real and imaginary part 
    index = @(j, k) (j - 1) * N + k - j - j * (j - 1) / 2;
    next_row = -1;
    for j = 1 : N - 1
        for k = j + 1 : N
            next_row = next_row + 2; % start with index = 1
            Ukk_Ujk = U(k, k) * U(j, k)';
            R(next_row + 0, index(j, k)) = - 2 * real(Ukk_Ujk);
            R(next_row + 1, index(j, k)) = - 2 * imag(Ukk_Ujk);
            for l = 1 : N
                if (l ~= j && l ~= k)
                    if (k < l)
                        Ukl = U(k, l);
                        index_kl = index(k, l);
                    else
                        Ukl = U(l, k)';
                        index_kl = index(l, k);
                    end
                    if (j < l)
                        Ujl = U(j, l)';
                        index_jl = index(j, l);
                    else
                        Ujl = U(l, j);
                        index_jl = index(l, j);
                    end
                    Ukl_Ujl = Ukl * Ujl;
                    UR = real(Ukl_Ujl);
                    UI = imag(Ukl_Ujl);
                    R(next_row + 0, index_kl) = UR * sign(l - k);
                    R(next_row + 1, index_kl) = UI * sign(l - k);
                    R(next_row + 0, index_jl) = UR * sign(j - l);
                    R(next_row + 1, index_jl) = UI * sign(j - l);
                end
            end
        end
    end

    % [U S V] = svd(R); 
    % svd(R)
    % sum_of_S = sum(S.^2)
    % trace_R = trace(R'*R)
    % return

    tau
    trivial_phases
    z
    switch METHOD
        case 'S'
            if VERBOSE_MODE
                disp(sprintf('Method: ''S'' with ''SV_TOLERANCE'' = %g', SV_TOLERANCE));
                disp(sprintf('Wait... Getting SVD of ''R'''));
            end
            NZSV_R = sum((svd(R) > SV_TOLERANCE)); % number of non-zero singular values of R
            d = tau - NZSV_R - trivial_phases - z;
        case 'R'
            if VERBOSE_MODE
                disp(sprintf('Method: ''R'''));
                disp(sprintf('Wait... Getting rank of ''R'''));
            end
            d = tau - rank(R) - trivial_phases - z;
        otherwise
            error('Not implemented!');
    end
end


function verbose=VERBOSE_MODE
% Toggle:
%   verbose = false  - to suppress all displays except for the defect value
%   verbose = true   - to display everything (helpful when dealing with really big matrices, to get a hint what's currently going on)
    verbose = false;
end

