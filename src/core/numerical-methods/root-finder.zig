const std = @import("std");

pub const RootFinderError = error{
    ToleranceBelowZero,
    MaxIterationsReached,
};

fn solver(
    f: anytype,
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64,
    tries: usize,
    max_iter: usize,
    tol: f64,
) RootFinderError!f64 {
    if (@abs(y1) <= tol) return x1;
    if (tries >= max_iter) return RootFinderError.MaxIterationsReached;

    var x2: f64 = undefined;
    if (@abs(y1) > @abs(y0)) {
        if (@abs(y1) <= std.math.floatEps(f64)) return x1;
        const ratio = y0 / y1;
        x2 = (-ratio * x1 + x0) / (1.0 - ratio);
    } else {
        if (@abs(y0) <= std.math.floatEps(f64)) return x0;
        const ratio = y1 / y0;
        x2 = (-ratio * x0 + x1) / (1.0 - ratio);
    }

    const y2 = f(x2);
    if (@abs(y2) <= tol) return x2;

    return solver(f, x1, y1, x2, y2, tries + 1, max_iter, tol);
}

pub fn secantMethod(f: anytype, x0: f64, tol: f64) RootFinderError!f64 {
    const max_iter: usize = 50;
    if (tol < 0.0) return RootFinderError.ToleranceBelowZero;

    const eps: f64 = 1e-4;
    const x1_shift = x0 * (1.0 + eps);
    const x1 = x1_shift + if (x1_shift >= 0.0) eps else -eps;

    const y0 = f(x0);
    if (@abs(y0) <= tol) return x0;

    const y1 = f(x1);
    return solver(f, x0, y0, x1, y1, 0, max_iter, tol);
}

test "secant method polynomial root" {
    const f = struct {
        fn call(x: f64) f64 {
            return x * x - 2.0 * x - 1.0;
        }
    };

    const root = try secantMethod(f.call, 3.0, 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), f.call(root), 1e-6);
}

test "secant method transcendental root" {
    const f = struct {
        fn call(x: f64) f64 {
            return std.math.exp(x) - std.math.sin(x);
        }
    };

    const root = try secantMethod(f.call, 3.0, 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), f.call(root), 1e-6);
}

test "secant method negative tolerance" {
    const f = struct {
        fn call(x: f64) f64 {
            return x;
        }
    };

    try std.testing.expectError(RootFinderError.ToleranceBelowZero, secantMethod(f.call, 3.0, -1e-6));
}

test "secant method fails to converge" {
    const f = struct {
        fn call(_: f64) f64 {
            return 1.0;
        }
    };

    try std.testing.expectError(RootFinderError.MaxIterationsReached, secantMethod(f.call, 3.0, 1e-15));
}
