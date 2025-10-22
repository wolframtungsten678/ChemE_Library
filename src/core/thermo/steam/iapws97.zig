const std = @import("std");

const units = @import("../../units.zig");
const constants = @import("../../constants.zig");
const steam_constants = @import("iapws97-constants.zig");
const root_finder = @import("../../numerical-methods/root-finder.zig");
const IjnRegionPoint = steam_constants.IjnRegionPoint;
const JnRegionPoint = steam_constants.JnRegionPoint;

const K = units.K;
const Pa = units.Pa;
const JPerKgK = units.JPerKgK;
const Pressure = units.Pressure;
const Temperature = units.Temperature;
const EnergyPerMass = units.EnergyPerMass;
const EnergyPerMassTemperature = units.EnergyPerMassTemperature;
const Velocity = units.Velocity;
const SpecificVolume = units.SpecificVolume;
const JPerKg = units.JPerKg;
const MPerSec = units.MPerSec;
const M3PerKg = units.M3PerKg;
const GAS_CONSTANT = constants.WATER_GAS_CONSTANT;
const SteamError = error{
    AboveCriticalTemperature,
    AboveCriticalPressure,
    BelowCriticalTemperature,
    TemperatureLow,
    TemperatureHigh,
    PressureLow,
    PressureHigh,
    FailedToConverge,
    InvalidPhaseFractions,
};

pub const SteamNonCriticalPhaseRegion = enum {
    Liquid,
    Vapor,
};

pub const LiquidVapor = struct {
    liquid_frac: f64,
    vapor_frac: f64,

    pub fn new(liquid_frac: f64, vapor_frac: f64) SteamError!LiquidVapor {
        if (liquid_frac < 0.0 or vapor_frac < 0.0) return SteamError.InvalidPhaseFractions;
        const sum = liquid_frac + vapor_frac;
        if (!std.math.approxEqAbs(f64, sum, 1.0, 1e-6)) return SteamError.InvalidPhaseFractions;
        return LiquidVapor{
            .liquid_frac = liquid_frac,
            .vapor_frac = vapor_frac,
        };
    }

    pub fn liquidFraction(self: LiquidVapor) f64 {
        return self.liquid_frac;
    }

    pub fn vaporFraction(self: LiquidVapor) f64 {
        return self.vapor_frac;
    }

    pub fn getLiquidFrac(self: LiquidVapor) f64 {
        return self.liquid_frac;
    }

    pub fn getVaporFrac(self: LiquidVapor) f64 {
        return self.vapor_frac;
    }
};

pub const CompositePhaseRegion = struct {
    liquid_vapor: LiquidVapor,

    pub fn liquidVapor(self: CompositePhaseRegion) LiquidVapor {
        return self.liquid_vapor;
    }
};

pub const PtPoint = struct {
    pressure: Pressure,
    temperature: Temperature,
};

pub const SatQuery = union(enum) {
    SatTQuery: struct {
        temperature: Temperature,
        phase_region: SteamNonCriticalPhaseRegion,
    },
    SatPQuery: struct {
        pressure: Pressure,
        phase_region: SteamNonCriticalPhaseRegion,
    },
};

pub const SteamQuery = union(enum) {
    Pt: PtPoint,
    Sat: SatQuery,
    EntropyP: struct {
        entropy: EnergyPerMassTemperature,
        pressure: Pressure,
    },
    EnthalpyP: struct {
        enthalpy: EnergyPerMass,
        pressure: Pressure,
    },
};

pub const PtvEntry = struct {
    temperature: Temperature,
    pressure: Pressure,
    phase_region: PhaseKind,
    internal_energy: EnergyPerMass,
    enthalpy: EnergyPerMass,
    entropy: EnergyPerMassTemperature,
    cv: EnergyPerMassTemperature,
    cp: EnergyPerMassTemperature,
    speed_of_sound: Velocity,
    specific_volume: SpecificVolume,
};

pub const PhaseKind = union(enum) {
    SupercriticalFluid: void,
    Gas: void,
    NonCritical: NonCriticalPhaseRegion,
    Composite: CompositePhaseRegion,
};

pub const NonCriticalPhaseRegion = enum {
    Liquid,
    Vapor,
};

const SpecificRegionPoint = struct {
    point: PtPoint,
    tau: f64,
    pi: f64,
    gamma: f64,
    gamma_pi: f64,
    gamma_pi_pi: f64,
    gamma_tau: f64,
    gamma_tau_tau: f64,
    gamma_pi_tau: f64,
};

var region3_current_point: PtPoint = .{
    .pressure = Pressure{ .pa = Pa.init(0.0) },
    .temperature = Temperature{ .k = K.init(0.0) },
};

fn pow(comptime T: type, base: T, exponent: f64) T {
    return std.math.pow(T, base, exponent);
}

fn powi(comptime T: type, base: T, exponent: i32) T {
    return std.math.pow(T, base, @as(T, @floatFromInt(exponent)));
}

fn getSatPressure(temperature: K) SteamError!Pa {
    const t = temperature.value;
    if (t > constants.WATER_CRITICAL_TEMPERATURE.value) return SteamError.AboveCriticalTemperature;

    const sat_temp_ratio = t / 1.0;
    const theta = sat_temp_ratio + (steam_constants.REGION_4[8].n /
        (sat_temp_ratio - steam_constants.REGION_4[9].n));
    const a = powi(f64, theta, 2) + steam_constants.REGION_4[0].n * theta +
        steam_constants.REGION_4[1].n;
    const b = steam_constants.REGION_4[2].n * powi(f64, theta, 2) +
        steam_constants.REGION_4[3].n * theta + steam_constants.REGION_4[4].n;
    const c = steam_constants.REGION_4[5].n * powi(f64, theta, 2) +
        steam_constants.REGION_4[6].n * theta + steam_constants.REGION_4[7].n;
    const inner = powi(f64, b, 2) - 4.0 * a * c;
    const pressure = powi(
        f64,
        (2.0 * c) / (-b + std.math.sqrt(inner)),
        4,
    ) * 1e6;
    return Pa.init(pressure);
}

fn getSatTemperature(pressure: Pa) SteamError!K {
    const p = pressure.value;
    if (p > constants.WATER_CRITICAL_PRESSURE.value) return SteamError.AboveCriticalPressure;

    const beta = std.math.pow(f64, p / 1e6, 0.25);
    const e = powi(f64, beta, 2) + steam_constants.REGION_4[2].n * beta +
        steam_constants.REGION_4[5].n;
    const f = steam_constants.REGION_4[0].n * powi(f64, beta, 2) +
        steam_constants.REGION_4[3].n * beta + steam_constants.REGION_4[6].n;
    const g = steam_constants.REGION_4[1].n * powi(f64, beta, 2) +
        steam_constants.REGION_4[4].n * beta + steam_constants.REGION_4[7].n;
    const d = (2.0 * g) / (-f - std.math.sqrt(powi(f64, f, 2) - 4.0 * e * g));
    const numerator = steam_constants.REGION_4[9].n + d -
        std.math.sqrt(powi(f64, steam_constants.REGION_4[9].n + d, 2) -
            4.0 * (steam_constants.REGION_4[8].n + steam_constants.REGION_4[9].n * d));
    const temperature = numerator / 2.0;
    return K.init(temperature);
}

fn getBoundary34Pressure(temperature: K) SteamError!Pa {
    const t = temperature.value;
    if (t < constants.WATER_CRITICAL_TEMPERATURE.value) return SteamError.BelowCriticalTemperature;

    const theta = t / 1.0;
    const pressure = (steam_constants.BOUNDARY_34[0].n +
        steam_constants.BOUNDARY_34[1].n * theta +
        steam_constants.BOUNDARY_34[2].n * powi(f64, theta, 2)) * 1e6;
    return Pa.init(pressure);
}

fn extractPressure(query: SteamQuery) ?Pressure {
    return switch (query) {
        .Pt => |pt| pt.pressure,
        .Sat => |sat| switch (sat) {
            .SatPQuery => |sat_p| sat_p.pressure,
            else => null,
        },
        .EntropyP => |data| data.pressure,
        .EnthalpyP => |data| data.pressure,
    };
}

fn extractTemperature(query: SteamQuery) ?Temperature {
    return switch (query) {
        .Pt => |pt| pt.temperature,
        .Sat => |sat| switch (sat) {
            .SatTQuery => |sat_t| sat_t.temperature,
            else => null,
        },
        else => null,
    };
}

fn checkIfOutOfRange(query: SteamQuery) SteamError!void {
    const opt_p = if (extractPressure(query)) |p| p.convertToSiUnit().value else null;
    const opt_t = if (extractTemperature(query)) |t| t.convertToSiUnit().value else null;

    if (opt_t) |t| {
        if (t < 273.15) return SteamError.TemperatureLow;
        if (t > 2000.0 + 273.15) return SteamError.TemperatureHigh;
    }

    if (opt_p) |p_val| {
        if (opt_t) |t_val| {
            if (p_val > 50e6 and t_val > 800.0 + 273.15) return SteamError.TemperatureHigh;
        }
        if (p_val < 0.0) return SteamError.PressureLow;
        if (p_val > 100e6) return SteamError.PressureHigh;
    }
}

const Iapws97Region = enum {
    Region1,
    Region2,
    Region3,
    Region4,
    Region5,
};

fn getRegionFromPtPoint(pt_point: PtPoint) SteamError!Iapws97Region {
    const t_si = pt_point.temperature.convertToSiUnit();
    const p = pt_point.pressure.convertToSiUnit().value;
    const t = t_si.value;

    const sat_pressure = getSatPressure(t_si) catch |err| switch (err) {
        SteamError.AboveCriticalTemperature => return Iapws97Region.Region5,
        else => return err,
    };

    const boundary_pressure = getBoundary34Pressure(t_si) catch |err| switch (err) {
        SteamError.BelowCriticalTemperature => null,
        else => return err,
    };

    if (t > 273.15 + 800.0) return Iapws97Region.Region5;
    if (t > 273.15 + 600.0) return Iapws97Region.Region2;

    if (p == sat_pressure.value) return Iapws97Region.Region4;
    if (p < sat_pressure.value) return Iapws97Region.Region2;
    if (boundary_pressure) |boundary| {
        if (p < boundary.value) return Iapws97Region.Region2;
        return Iapws97Region.Region3;
    }
    return Iapws97Region.Region1;
}

fn getRegionFromSatQuery(sat_query: SatQuery) SteamError!struct { pt: PtPoint, region: Iapws97Region } {
    const phase_region: SteamNonCriticalPhaseRegion = switch (sat_query) {
        .SatTQuery => |sat| sat.phase_region,
        .SatPQuery => |sat| sat.phase_region,
    };

    return switch (sat_query) {
        .SatTQuery => |sat| blk: {
            const p = try getSatPressure(sat.temperature.convertToSiUnit());
            break :blk .{
                .pt = PtPoint{
                    .pressure = Pressure{ .pa = p },
                    .temperature = sat.temperature,
                },
                .region = if (phase_region == .Liquid) Iapws97Region.Region1 else Iapws97Region.Region2,
            };
        },
        .SatPQuery => |sat| blk: {
            const t = try getSatTemperature(sat.pressure.convertToSiUnit());
            break :blk .{
                .pt = PtPoint{
                    .pressure = sat.pressure,
                    .temperature = Temperature{ .k = t },
                },
                .region = if (phase_region == .Liquid) Iapws97Region.Region1 else Iapws97Region.Region2,
            };
        },
    };
}

fn createEntryFromRegionPoint(
    specific_region_point: SpecificRegionPoint,
    phase_region: PhaseKind,
) PtvEntry {
    const temperature = specific_region_point.point.temperature.convertToSiUnit().value;
    const pressure = specific_region_point.point.pressure.convertToSiUnit().value;
    const pi = specific_region_point.pi;
    const tau = specific_region_point.tau;
    const gamma = specific_region_point.gamma;
    const gamma_pi = specific_region_point.gamma_pi;
    const gamma_pi_pi = specific_region_point.gamma_pi_pi;
    const gamma_tau = specific_region_point.gamma_tau;
    const gamma_tau_tau = specific_region_point.gamma_tau_tau;
    const gamma_pi_tau = specific_region_point.gamma_pi_tau;

    const internal_energy = GAS_CONSTANT.value * temperature * (tau * gamma_tau - pi * gamma_pi);
    const enthalpy = GAS_CONSTANT.value * temperature * tau * gamma_tau;
    const entropy = GAS_CONSTANT.value * (tau * gamma_tau - gamma);
    const cv = GAS_CONSTANT.value *
        (-std.math.pow(f64, tau, 2.0) * gamma_tau_tau +
            std.math.pow(f64, gamma_pi - tau * gamma_pi_tau, 2.0) / gamma_pi_pi);
    const cp = GAS_CONSTANT.value * -std.math.pow(f64, tau, 2.0) * gamma_tau_tau;
    const speed_of_sound = std.math.sqrt(GAS_CONSTANT.value * temperature *
        (std.math.pow(f64, gamma_pi, 2.0) /
            ((std.math.pow(f64, gamma_pi - tau * gamma_pi_tau, 2.0) /
                (std.math.pow(f64, tau, 2.0) * gamma_tau_tau)) - gamma_pi_pi)));
    const specific_volume = pi * (gamma_pi * GAS_CONSTANT.value * temperature) / pressure;

    return PtvEntry{
        .temperature = Temperature{ .k = K.init(temperature) },
        .pressure = Pressure{ .pa = Pa.init(pressure) },
        .phase_region = phase_region,
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(internal_energy) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(enthalpy) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(entropy) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(cv) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(cp) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(speed_of_sound) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(specific_volume) },
    };
}

fn gibbsMethod(point: PtPoint) PtvEntry {
    const pressure = point.pressure.convertToSiUnit().value;
    const temperature = point.temperature.convertToSiUnit().value;
    const pi = pressure / 16.53e6;
    const tau = 1386.0 / temperature;

    var gamma: f64 = 0.0;
    var gamma_pi: f64 = 0.0;
    var gamma_pi_pi: f64 = 0.0;
    var gamma_tau: f64 = 0.0;
    var gamma_tau_tau: f64 = 0.0;
    var gamma_pi_tau: f64 = 0.0;

    for (steam_constants.REGION_1_AND_4) |region_point| {
        const n = region_point.n;
        const i = region_point.i;
        const j = region_point.j;
        const pow_pi = pow(f64, 7.1 - pi, i);
        const pow_tau = pow(f64, tau - 1.222, j);
        gamma += n * pow_pi * pow_tau;
        gamma_pi += -n * i * pow(f64, 7.1 - pi, i - 1.0) * pow_tau;
        gamma_pi_pi += n * i * (i - 1.0) * pow(f64, 7.1 - pi, i - 2.0) * pow_tau;
        gamma_tau += n * j * pow_pi * pow(f64, tau - 1.222, j - 1.0);
        gamma_tau_tau += n * j * (j - 1.0) * pow_pi * pow(f64, tau - 1.222, j - 2.0);
        gamma_pi_tau += -n * i * j * pow(f64, 7.1 - pi, i - 1.0) * pow(f64, tau - 1.222, j - 1.0);
    }

    const specific_region_point = SpecificRegionPoint{
        .point = point,
        .tau = tau,
        .pi = pi,
        .gamma = gamma,
        .gamma_pi = gamma_pi,
        .gamma_pi_pi = gamma_pi_pi,
        .gamma_tau = gamma_tau,
        .gamma_tau_tau = gamma_tau_tau,
        .gamma_pi_tau = gamma_pi_tau,
    };

    return createEntryFromRegionPoint(specific_region_point, PhaseKind{ .NonCritical = .Liquid });
}

fn vaporMethod(
    tau: f64,
    tau_shift: f64,
    point: PtPoint,
    ideal_points: []const JnRegionPoint,
    residual_points: []const IjnRegionPoint,
) PtvEntry {
    const pressure = point.pressure.convertToSiUnit().value;
    const temperature = point.temperature.convertToSiUnit().value;
    const pi = pressure / 1.0e6;

    var gamma = @log(pi);
    var gamma_pi = 1.0 / pi;
    var gamma_pi_pi = -1.0 / std.math.pow(f64, pi, 2.0);
    var gamma_tau: f64 = 0.0;
    var gamma_tau_tau: f64 = 0.0;
    var gamma_pi_tau: f64 = 0.0;

    for (ideal_points) |region_point| {
        const n = region_point.n;
        const j = region_point.j;
        gamma += n * pow(f64, tau, j);
        gamma_tau += n * j * pow(f64, tau, j - 1.0);
        gamma_tau_tau += n * j * (j - 1.0) * pow(f64, tau, j - 2.0);
    }

    for (residual_points) |region_point| {
        const n = region_point.n;
        const i = region_point.i;
        const j = region_point.j;
        const pow_pi = pow(f64, pi, i);
        const pow_tau = pow(f64, tau - tau_shift, j);
        gamma += n * pow_pi * pow_tau;
        gamma_pi += n * i * pow(f64, pi, i - 1.0) * pow_tau;
        gamma_pi_pi += n * i * (i - 1.0) * pow(f64, pi, i - 2.0) * pow_tau;
        gamma_tau += n * pow_pi * j * pow(f64, tau - tau_shift, j - 1.0);
        gamma_tau_tau += n * pow_pi * j * (j - 1.0) * pow(f64, tau - tau_shift, j - 2.0);
        gamma_pi_tau += n * i * pow(f64, pi, i - 1.0) * j * pow(f64, tau - tau_shift, j - 1.0);
    }

    const phase_info = blk: {
        const above_temp = temperature > constants.WATER_CRITICAL_TEMPERATURE.value;
        const above_pressure = pressure > constants.WATER_CRITICAL_PRESSURE.value;
        if (above_temp and above_pressure) {
            break :blk PhaseKind{ .SupercriticalFluid = {} };
        } else if (above_temp and !above_pressure) {
            break :blk PhaseKind{ .Gas = {} };
        } else {
            break :blk PhaseKind{ .NonCritical = .Vapor };
        }
    };

    const specific_region_point = SpecificRegionPoint{
        .point = point,
        .tau = tau,
        .pi = pi,
        .gamma = gamma,
        .gamma_pi = gamma_pi,
        .gamma_pi_pi = gamma_pi_pi,
        .gamma_tau = gamma_tau,
        .gamma_tau_tau = gamma_tau_tau,
        .gamma_pi_tau = gamma_pi_tau,
    };

    return createEntryFromRegionPoint(specific_region_point, phase_info);
}

fn region3BySpecificVolume(pt_point: PtPoint, specific_volume: f64) PtvEntry {
    const density = 1.0 / specific_volume;
    const n1 = steam_constants.REGION_3_N1.n;
    const delta = density / 322.0;
    const temperature = pt_point.temperature.convertToSiUnit().value;
    const tau = 647.096 / temperature;

    var phi = n1 * @log(delta);
    var phi_delta = n1 / delta;
    var phi_delta_delta = -n1 / std.math.pow(f64, delta, 2.0);
    var phi_tau: f64 = 0.0;
    var phi_tau_tau: f64 = 0.0;
    var phi_delta_tau: f64 = 0.0;

    for (steam_constants.REGION_3) |region_point| {
        const n = region_point.n;
        const i = region_point.i;
        const j = region_point.j;
        const pow_delta = pow(f64, delta, i);
        const pow_tau = pow(f64, tau, j);
        phi += n * pow_delta * pow_tau;
        phi_delta += n * i * pow(f64, delta, i - 1.0) * pow_tau;
        phi_delta_delta += n * i * (i - 1.0) * pow(f64, delta, i - 2.0) * pow_tau;
        phi_tau += n * pow_delta * j * pow(f64, tau, j - 1.0);
        phi_tau_tau += n * pow_delta * j * (j - 1.0) * pow(f64, tau, j - 2.0);
        phi_delta_tau += n * i * pow(f64, delta, i - 1.0) * j * pow(f64, tau, j - 1.0);
    }

    const pressure = phi_delta * delta * density * GAS_CONSTANT.value * temperature;
    const internal_energy = tau * phi_tau * GAS_CONSTANT.value * temperature;
    const enthalpy = (tau * phi_tau + delta * phi_delta) * GAS_CONSTANT.value * temperature;
    const entropy = (tau * phi_tau - phi) * GAS_CONSTANT.value;
    const cv = -std.math.pow(f64, tau, 2.0) * phi_tau_tau * GAS_CONSTANT.value;
    const cp = (-std.math.pow(f64, tau, 2.0) * phi_tau_tau +
        std.math.pow(f64, delta * phi_delta - delta * tau * phi_delta_tau, 2.0) /
            (2.0 * delta * phi_delta + std.math.pow(f64, delta, 2.0) * phi_delta_delta)) *
        GAS_CONSTANT.value;

    const speed_of_sound = std.math.sqrt(
        (2.0 * delta * phi_delta + std.math.pow(f64, delta, 2.0) * phi_delta_delta -
            std.math.pow(f64, delta * phi_delta - delta * tau * phi_delta_tau, 2.0) /
                (std.math.pow(f64, tau, 2.0) * phi_tau_tau)) * GAS_CONSTANT.value * temperature,
    );

    return PtvEntry{
        .temperature = Temperature{ .k = K.init(temperature) },
        .pressure = Pressure{ .pa = Pa.init(pressure) },
        .phase_region = PhaseKind{ .SupercriticalFluid = {} },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(internal_energy) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(enthalpy) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(entropy) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(cv) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(cp) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(speed_of_sound) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(specific_volume) },
    };
}

fn region3Residual(x: f64) f64 {
    const entry = region3BySpecificVolume(region3_current_point, x);
    return entry.pressure.convertToSiUnit().value -
        region3_current_point.pressure.convertToSiUnit().value;
}

fn region3Method(point: PtPoint) SteamError!PtvEntry {
    region3_current_point = point;
    const specific_volume = root_finder.secantMethod(region3Residual, 1.0 / 500.0, 1e-4) catch {
        return SteamError.FailedToConverge;
    };
    return region3BySpecificVolume(point, specific_volume);
}

fn getEntryFromPtPoint(point: PtPoint, region: Iapws97Region) SteamError!PtvEntry {
    const temperature = point.temperature.convertToSiUnit().value;
    return switch (region) {
        .Region1, .Region4 => gibbsMethod(point),
        .Region2 => vaporMethod(
            540.0 / temperature,
            0.5,
            point,
            steam_constants.REGION_2_IDEAL,
            steam_constants.REGION_2_RESIDUAL,
        ),
        .Region3 => try region3Method(point),
        .Region5 => vaporMethod(
            1000.0 / temperature,
            0.0,
            point,
            steam_constants.REGION_5_IDEAL,
            steam_constants.REGION_5_RESIDUAL,
        ),
    };
}

fn interpolateEntry(
    liquid_entry: PtvEntry,
    vapor_entry: PtvEntry,
    liq_frac: f64,
) SteamError!PtvEntry {
    const vap_frac = 1.0 - liq_frac;
    const interpolate = struct {
        fn apply(
            liquid: PtvEntry,
            vapor: PtvEntry,
            liq_fraction: f64,
            vap_fraction: f64,
            comptime accessor: fn (PtvEntry) f64,
        ) f64 {
            return accessor(liquid) * liq_fraction + accessor(vapor) * vap_fraction;
        }
    };

    const phase_info = LiquidVapor.new(liq_frac, vap_frac) catch |err| {
        return err;
    };

    const temperature = interpolate.apply(
        liquid_entry,
        vapor_entry,
        liq_frac,
        vap_frac,
        struct {
            fn value(entry: PtvEntry) f64 {
                return entry.temperature.convertToSiUnit().value;
            }
        }.value,
    );

    const pressure = interpolate.apply(
        liquid_entry,
        vapor_entry,
        liq_frac,
        vap_frac,
        struct {
            fn value(entry: PtvEntry) f64 {
                return entry.pressure.convertToSiUnit().value;
            }
        }.value,
    );

    const internal_energy = interpolate.apply(
        liquid_entry,
        vapor_entry,
        liq_frac,
        vap_frac,
        struct {
            fn value(entry: PtvEntry) f64 {
                return entry.internal_energy.convertToSiUnit().value;
            }
        }.value,
    );

    const enthalpy = interpolate.apply(
        liquid_entry,
        vapor_entry,
        liq_frac,
        vap_frac,
        struct {
            fn value(entry: PtvEntry) f64 {
                return entry.enthalpy.convertToSiUnit().value;
            }
        }.value,
    );

    const entropy = interpolate.apply(
        liquid_entry,
        vapor_entry,
        liq_frac,
        vap_frac,
        struct {
            fn value(entry: PtvEntry) f64 {
                return entry.entropy.convertToSiUnit().value;
            }
        }.value,
    );

    const cv = interpolate.apply(
        liquid_entry,
        vapor_entry,
        liq_frac,
        vap_frac,
        struct {
            fn value(entry: PtvEntry) f64 {
                return entry.cv.convertToSiUnit().value;
            }
        }.value,
    );

    const cp = interpolate.apply(
        liquid_entry,
        vapor_entry,
        liq_frac,
        vap_frac,
        struct {
            fn value(entry: PtvEntry) f64 {
                return entry.cp.convertToSiUnit().value;
            }
        }.value,
    );

    const speed_of_sound = interpolate.apply(
        liquid_entry,
        vapor_entry,
        liq_frac,
        vap_frac,
        struct {
            fn value(entry: PtvEntry) f64 {
                return entry.speed_of_sound.convertToSiUnit().value;
            }
        }.value,
    );

    const specific_volume_recip = interpolate.apply(
        liquid_entry,
        vapor_entry,
        liq_frac,
        vap_frac,
        struct {
            fn value(entry: PtvEntry) f64 {
                return 1.0 / entry.specific_volume.convertToSiUnit().value;
            }
        }.value,
    );

    return PtvEntry{
        .temperature = Temperature{ .k = K.init(temperature) },
        .pressure = Pressure{ .pa = Pa.init(pressure) },
        .phase_region = PhaseKind{
            .Composite = CompositePhaseRegion{ .liquid_vapor = phase_info },
        },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(internal_energy) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(enthalpy) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(entropy) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(cv) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(cp) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(speed_of_sound) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(1.0 / specific_volume_recip) },
    };
}

fn iteratePtEntrySolution(
    pressure: Pa,
    target_value: f64,
    comptime getPropValue: fn (PtvEntry) f64,
) SteamError!PtvEntry {
    const liquid_query = SteamQuery{
        .Sat = SatQuery{
            .SatPQuery = .{
                .pressure = Pressure{ .pa = pressure },
                .phase_region = SteamNonCriticalPhaseRegion.Liquid,
            },
        },
    };

    const vapor_query = SteamQuery{
        .Sat = SatQuery{
            .SatPQuery = .{
                .pressure = Pressure{ .pa = pressure },
                .phase_region = SteamNonCriticalPhaseRegion.Vapor,
            },
        },
    };

    const liquid_entry = try getSteamTableEntry(liquid_query);
    const vapor_entry = try getSteamTableEntry(vapor_query);

    const liquid_value = getPropValue(liquid_entry);
    const vapor_value = getPropValue(vapor_entry);
    if (liquid_value <= target_value and vapor_value >= target_value) {
        const liq_frac = (vapor_value - target_value) /
            (vapor_value - liquid_value);
        return interpolateEntry(liquid_entry, vapor_entry, liq_frac);
    }

    const ctx = struct {
        pressure: Pa,
        target: f64,

        fn eval(self: @This(), temperature: f64) f64 {
            const query = SteamQuery{
                .Pt = PtPoint{
                    .pressure = Pressure{ .pa = self.pressure },
                    .temperature = Temperature{ .k = K.init(temperature) },
                },
            };
            const result = getSteamTableEntry(query) catch return std.math.nan(f64);
            return getPropValue(result) - self.target;
        }
    }{ .pressure = pressure, .target = target_value };

    const tol = 1e-5;
    if (tol < 0.0) return SteamError.FailedToConverge;
    const max_iter: usize = 50;
    const eps: f64 = 1e-4;

    var x0: f64 = 310.0;
    var y0 = ctx.eval(x0);
    if (!std.math.isFinite(y0)) return SteamError.FailedToConverge;
    if (@abs(y0) <= tol) {
        return getSteamTableEntry(SteamQuery{
            .Pt = PtPoint{
                .pressure = Pressure{ .pa = pressure },
                .temperature = Temperature{ .k = K.init(x0) },
            },
        });
    }

    const x1_shift: f64 = x0 * (1.0 + eps);
    var x1: f64 = x1_shift + if (x1_shift >= 0.0) eps else -eps;
    var y1 = ctx.eval(x1);
    if (!std.math.isFinite(y1)) return SteamError.FailedToConverge;

    var tries: usize = 0;
    while (tries < max_iter) : (tries += 1) {
        if (@abs(y1) <= tol) {
            const final_query = SteamQuery{
                .Pt = PtPoint{
                    .pressure = Pressure{ .pa = pressure },
                    .temperature = Temperature{ .k = K.init(x1) },
                },
            };
            return getSteamTableEntry(final_query);
        }

        var x2: f64 = undefined;
        if (@abs(y1) > @abs(y0)) {
            if (@abs(y1) <= std.math.floatEps(f64)) break;
            const ratio = y0 / y1;
            x2 = (-ratio * x1 + x0) / (1.0 - ratio);
        } else {
            if (@abs(y0) <= std.math.floatEps(f64)) break;
            const ratio = y1 / y0;
            x2 = (-ratio * x0 + x1) / (1.0 - ratio);
        }

        const y2 = ctx.eval(x2);
        if (!std.math.isFinite(y2)) return SteamError.FailedToConverge;
        if (@abs(y2) <= tol) {
            const final_query = SteamQuery{
                .Pt = PtPoint{
                    .pressure = Pressure{ .pa = pressure },
                    .temperature = Temperature{ .k = K.init(x2) },
                },
            };
            return getSteamTableEntry(final_query);
        }

        x0 = x1;
        y0 = y1;
        x1 = x2;
        y1 = y2;
    }

    return SteamError.FailedToConverge;
}

pub fn getSteamTableEntry(query: SteamQuery) SteamError!PtvEntry {
    try checkIfOutOfRange(query);

    return switch (query) {
        .Pt => |point| blk: {
            const region = try getRegionFromPtPoint(point);
            std.debug.print("getRegionFromPtPoint {s}", .{@tagName(region)});
            break :blk try getEntryFromPtPoint(point, region);
        },
        .Sat => |sat_query| blk: {
            const result = try getRegionFromSatQuery(sat_query);
            std.debug.print("getRegionFromSatQuery {s}", .{@tagName(result.region)});
            break :blk try getEntryFromPtPoint(result.pt, result.region);
        },
        .EntropyP => |data| try iteratePtEntrySolution(
            data.pressure.convertToSiUnit(),
            data.entropy.convertToSiUnit().value,
            struct {
                fn get(entry: PtvEntry) f64 {
                    return entry.entropy.convertToSiUnit().value;
                }
            }.get,
        ),
        .EnthalpyP => |data| try iteratePtEntrySolution(
            data.pressure.convertToSiUnit(),
            data.enthalpy.convertToSiUnit().value,
            struct {
                fn get(entry: PtvEntry) f64 {
                    return entry.enthalpy.convertToSiUnit().value;
                }
            }.get,
        ),
    };
}

test "Steam Test Region 1" {
    const query = SteamQuery{ .Pt = PtPoint{
        .temperature = Temperature{ .k = K.init(750) },
        .pressure = Pressure{ .pa = Pa.init(78.309563916917e3) },
    } };

    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(750) },
        .pressure = Pressure{ .pa = Pa.init(78.309563916917e3) },
        .phase_region = PhaseKind{ .SupercriticalFluid = {} },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(2102.069317626429e3) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(2258.688445460262e3) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(4.469719056217e3) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.71701677121e3) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(6.341653594791e3) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(760.696040876798) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(@as(f64, 1) / @as(f64, 500)) },
    };

    const actual = try getSteamTableEntry(query);

    //try std.testing.expectApproxEqAbs(expected.Pt.pressure, actual.Pt.pressure, 1e-9);
    try std.testing.expectApproxEqAbs(expected.pressure.getValue(), actual.pressure.getValue(), 1e-9);
}
