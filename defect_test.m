% defect "unit test"
%
% 20170103
% W. Bruzda, name[at]uj.edu.pl : name = w.bruzda
%
% http://chaos.if.uj.edu.pl/~karol/hadamard/
% https://github.com/matrix-toolbox/
%
% >> version % 9.1.0.441655 (R2016b)
% >> result = runtests('defect_test.m');
% >> rt = table(result) 

function tests = defect_test
    tests = functiontests(localfunctions);
end

function setup(testCase)
    addpath ../chm; % include Fourier...
end

function testUnitaryDefect(testCase)

    prime_number = [ 2 3 5 7 ];

    for p = prime_number
        F = kron(fourier(p, 'classic'), fourier(p, 'classic')); % tensor product of 'classic' Fouriers

        expected_value = (p - 1) * (p - 1) * (p + 1);

        d = [ ...
            expected_value, ...
            get_defect(F)
        ];
        disp(sprintf('Expected value: %d', expected_value));
        disp('------------------------------------------');

        if ~all(d == d(1))
            error('Result does not match expected value!');
        end
    end
end

function testRestrictedDefect(testCase)

    dimension = ...
        [ 2 3  4  5   6   7   8   9  10  11   12  13   14 ];

    expected_value = ... % of the restricted defect for hermitian Fourier matrix
        [ 0 4 21 36 112 120 273 352 576 540 1237 924 1632 ];

    selected_dimension = [ 2 4 5 8 10 ];

    for N = selected_dimension
        F = fourier(N, 'hermitian');
        d = [
            expected_value(N - 1), ...
            get_defect(F);
        ];
        disp(sprintf('Expected value: %d', expected_value(N - 1)));
        disp('------------------------------------------'); 

        if ~all(d == d(1))
            error('Result does not match expected value in dimension: %n', N);
        end
    end

end

function result = get_defect(U)

    result = {};

    result{1} = defect(U, 'R');                    % default
    result{2} = defect(U, 'R', 1e-11);             % default with custom SV_TOLERANCE
    result{3} = defect(U, 'R', 'wrong_tolerance'); % wrong SV_TOLERANCE
    result{4} = defect(U);                         % default
    result{5} = defect(U, 1e-11);                  % default - SV_TOLERANCE should be ignored
    result{6} = defect(U, 'S', 1e-11);             % SVD with custom SV_TOLERANCE
    result{7} = defect(U, 'S');                    % SVD with default SV_TOLERANCE

    % FIX it using verifyError() or by implementing 'T' method in defect_h.m :)

    % result{8} = defect(U, 'T');                  % tangent...
    % result{9} = defect(U, 'T', 1e-11);           % tangent... SV_TOLERANCE should be ignored

    result{10} = defect(U, 'X');                   % wrong method - should switch to default
    result{11} = defect(U, 'X', 1e-11);            % wrong method - should switch to default and ignore SV_TOLERANCE

    result = cell2mat(result);
    disp('------------------------------------------');
    disp('All non-negative integers should be equal!');
    disp(result);
    disp('------------------------------------------');

end

