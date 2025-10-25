const std = @import("std");

pub const RootFinderError = error{
    ToleranceBelowZero,
    MaxIterationsReached,
};

pub fn Closure(comptime T: type) type {
    return struct {
        ctx: T,
        func: *const fn (*T, f64) f64,

        pub fn call(self: *@This(), x: f64) f64 {
            return self.func(&self.ctx, x);
        }
    };
}

fn solver(
    comptime T: type,
    c: *Closure(T),
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

    const y2 = c.call(x2);
    if (@abs(y2) <= tol) return x2;

    return solver(T, c, x1, y1, x2, y2, tries + 1, max_iter, tol);
}

pub fn secantMethod(comptime T: type, c: *Closure(T), x0: f64, tol: f64) RootFinderError!f64 {
    const max_iter: usize = 50;
    if (tol < 0.0) return RootFinderError.ToleranceBelowZero;

    const eps: f64 = 1e-4;
    const x1_shift = x0 * (1.0 + eps);
    const x1 = x1_shift + if (x1_shift >= 0.0) eps else -eps;

    const y0 = c.call(x0);
    if (@abs(y0) <= tol) return x0;

    const y1 = c.call(x1);
    return solver(T, c, x0, y0, x1, y1, 0, max_iter, tol);
}

fn polynomialRoot(_: *void, x: f64) f64 {
    return x * x - 2.0 * x - 1.0;
}

test "secant method polynomial root" {
    var c = Closure(void){
        .ctx = {},
        .func = polynomialRoot,
    };

    const root = try secantMethod(void, &c, 3.0, 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), c.call(root), 1e-6);
}

fn transcendentalRoot(_: *void, x: f64) f64 {
    return std.math.exp(x) - std.math.sin(x);
}

test "secant method transcendental root" {
    var c = Closure(void){
        .ctx = {},
        .func = transcendentalRoot,
    };

    const root = try secantMethod(void, &c, 3.0, 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), c.call(root), 1e-6);
}

fn returnOne(_: *void, _: f64) f64 {
    return @as(f64, 1);
}

test "secant method negative tolerance" {
    var c = Closure(void){
        .ctx = {},
        .func = returnOne,
    };

    try std.testing.expectError(RootFinderError.ToleranceBelowZero, secantMethod(void, &c, 3.0, -1e-6));
}

test "secant method fails to converge" {
    var c = Closure(void){
        .ctx = {},
        .func = returnOne,
    };

    try std.testing.expectError(RootFinderError.MaxIterationsReached, secantMethod(void, &c, 3.0, 1e-15));
}
