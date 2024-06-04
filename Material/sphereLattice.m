function points = sphereLattice(dim, N)
% sphereLattice - Generate approximately equidistant points on a hpyersphere.
% 
% Arguments: 
% - dim:   (integer, scalar) number of dimensions of the hypersphere whose surface
%          will be of dimension dim-1.
% - N:     (integer, scalar) number of points to distribute on the sphere's surface.
% 
% Return value: 
% - points:  (Nxdim real) N points in "dim" dimensions at unit distnace from 0.
% 
% All functions were translated directly from the following Python toolbox authored
% by Erik Brinkman: https://github.com/erikbrinkman/fibonacci_lattice
% 
% Details can be found on
% https://stackoverflow.com/questions/57123194/how-to-distribute-points-evenly-on-the-surface-of-hyperspheres-in-higher-dimensi/59279721#59279721
% 
% 2024 - Gatien Clement, Institut Langevin, Paris, France

    if dim < 2
        error('Dimension must be greater than one, you provided %d.', dim);
    elseif N < 1
        error('You must request at least one point, you requested %d.', N);
    end
    cube = cube_lattice(dim - 1, N);
    points = ones(N, dim);
    points(:, 1) = points(:, 1) .* sin(2 * pi * cube(:, 1));
    points(:, 2) = points(:, 2) .* cos(2 * pi * cube(:, 1));
    pdegs = zeros(N, 1);
    HLP = log(pi) / 2;
    for d = 2:dim-1
        targets = cube(:, d);
        mult = exp(gammaln((d + 1) / 2) - gammaln(d / 2) - HLP);
        [~, order] = sort(targets);
        pdegs = inv_int_sin_ms(pdegs, targets(order) / mult, d - 1, 0, pi, 1e-10);
        points(order, 1:d) = points(order, 1:d) .* sin(pdegs);
        points(order, d + 1) = points(order, d + 1) .* cos(pdegs);
    end

end

function res = n_primes(n)
    % Return the n first prime numbers.
    if n > 100
        error('Only dimensions <101 are supported.'); 
    end
    p = primes(100);
    res = p(1:n);
end

function points = cube_lattice(dim, num_points)
    % Generate num_points points over the dim dimensional cube
    if dim < 1
    error('dimension must be greater than zero: %d', dim);
    elseif num_points < 1
    error('must request at least one point: %d', num_points);
    end
    rest = sqrt(n_primes(dim - 1));
    mults = [1 / num_points, rest];
    points = mod((mults .* (0:num_points - 1)'), 1.0);
end

function result = int_sin_m(x, m)
    % Computes the integral of sin^m(t) dt from 0 to x recursively
    cosx = cos(x);
    sinx = sin(x);
    sinx2 = sinx.^2;
    start = mod(m, 2);

    if start == 0
        result = x;
        sinxp = sinx;
    else
        result = 1 - cosx;
        sinxp = sinx2;
    end

    for p = start + 1:2:m-1
        result = p / (p + 1) * result - cosx * sinxp / (p + 1);
        sinxp = sinxp * sinx2;
    end
end

function mid = inv_int_sin_m(target, m, lower, upper, atol)
    % Returns func inverse of mult * integral of sin(x) ** m
    mid = (lower + upper) / 2;
    approx = int_sin_m(mid, m);
    while abs(approx - target) > atol
        if approx > target
            upper = mid;
        else
            lower = mid;
        end
        mid = (upper + lower) / 2;
        approx = int_sin_m(mid, m);
    end
end

function results = inv_int_sin_ms(results, targets, m, lower, upper, atol)
    % Returns func inverse of mult * integral of sin(x) ** m
    % inverse is accurate to an absolute tolerance of atol, and
    % must be monotonically increasing over the interval lower
    % to upper
    
    if isempty(targets)
        return;
    elseif isscalar(targets)
        results(1) = inv_int_sin_m(targets(1), m, lower, upper, atol);
    else
        mid = (lower + upper) / 2;
        approx = int_sin_m(mid, m);
        ind = find(targets > approx, 1);
        if isempty(ind)
            ind = numel(targets) + 1;
        end
        results(1:ind-1) = inv_int_sin_ms(results(1:ind-1), targets(1:ind-1), m, lower, mid, atol);
        results(ind:end) = inv_int_sin_ms(results(ind:end), targets(ind:end), m, mid, upper, atol);
    end
end
