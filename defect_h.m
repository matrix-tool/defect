% 20170209
% name[at]alumni.uj.edu.pl, name = w.bruzda
% https://github.com/matrix-tool/defect

% [Restricted] defect of a unitary matrix U subjected to additional constraints.

% Since "defect_h.m" is "internal" function for "defect.m" its arguments are assumed valid!

function d = defect_h(U, METHOD, SV_TOLERANCE)

    N = size(U, 1);
    disp('Wait... Preparing matrix ''R''...');
    z = 0; % number of zeros above the diagonal in U
    for j = 1 : N - 1
        for k = j + 1 : N
            z = z + (abs(U(j, k)) < 1e-12);
        end
    end
    % z = MUBs * N * (N - 1) / 2 for MUBs

    tau = N * (N - 1) / 2;
    trivial_phases = N - 1; % see equivalence relation...
    R = zeros(tau * 2, tau); % system of linear equations matrix split onto real and imaginary part 
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

    switch METHOD
        case 'S'
            disp(sprintf('Method: ''S'' with ''SV_TOLERANCE'' = %g', SV_TOLERANCE));
            disp(sprintf('Wait... Getting SVD of ''R'''));
            NZSV_R = sum((svd(R) > SV_TOLERANCE)); % number of non-zero singular values of R
            d = tau - NZSV_R - trivial_phases - z;
        case 'R'
            disp(sprintf('Method: ''R'''));
            disp(sprintf('Wait... Getting rank of ''R'''));
            d = tau - rank(R) - trivial_phases - z;
        otherwise
            error('Not implemented!');
    end
end

