% 20170304
% [dephased] defect of a square matrix

function d = defect(U, METHOD, SV_TOLERANCE)

    if ~exist('U') || ~ismatrix(U) || size(U, 1) ~= size(U, 2)
        error('No square matrix provided!');
    end

    N = size(U, 1);
    if N > 128
        warning(sprintf('Matrix size is %d. Higher orders take longer...', N));
    end

    try
        matrix_type = get_matrix_type(U);
    catch
        error('Unsupported matrix type!');
        return
    end

    % argument validator
    POSSIBLE_METHOD = ['R', 'S', 'T'];
    if ~exist('METHOD') || ~any(ismember(POSSIBLE_METHOD, METHOD))
        disp('No valid ''METHOD'' provided when calling ''defect(MATRIX [, METHOD [, SV_TOLERANCE]])''.');
        disp('Possible values:');
        disp('METHOD: ''R'' - rank of matrix ''R''.');
        disp('METHOD: ''S'' - number of non-zero singular values of ''R''.');
        disp('METHOD: ''T'' - dimension of tangent space...');
        disp('Attempt to calculate defect using the rank of ''R'' matrix.');
        METHOD = 'R';
    end
    if ~exist('SV_TOLERANCE') || ~isnumeric(SV_TOLERANCE)
        SV_TOLERANCE = 1e-13;
    end

    switch matrix_type
        case 'unitary'
            disp('Input matrix is (only) unitary.');
            disp('Continue with STANDARD defect procedure.');
            d = defect_u(U, METHOD, SV_TOLERANCE);
        case 'unitary_hermitian'
            disp('Input matrix is unitary and hermitian.');
            disp('Continue with RESTRICTED defect procedure.');
            d = defect_h(U, METHOD, SV_TOLERANCE);
        otherwise
            error('Not implemented...');
    end
end

function matrix_type = get_matrix_type(U)
    if is_unitary(U)
        matrix_type = 'unitary';
        if is_hermitian(U) % in-built "ishermitian" function has to small tolerance...
            matrix_type = 'unitary_hermitian';
        end
    end
end

function b = is_hermitian(U)
    b = false;
    T = abs(U - U');
    if sum(T(:)) < 1e-10 % note - the tolerance is relatively large...
        b = true;
    end
end

function b = is_unitary(U)
    b = false;
    T = U * U';
    T = T / T(1, 1); % U = unitary => T ~ I
    if norm(T - eye(size(T, 1)), 'fro') < 1e-10 % note - the tolerance is relatively large...
        b = true;
    end
end

