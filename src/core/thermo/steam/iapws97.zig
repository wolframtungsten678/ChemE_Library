const std = @import("std");

const units = @import("../../units.zig");
const constants = @import("../../constants.zig");
const steam_constants = @import("iapws97-constants.zig");
const root_finder = @import("../../numerical-methods/root-finder.zig");
const IjnRegionPoint = steam_constants.IjnRegionPoint;
const JnRegionPoint = steam_constants.JnRegionPoint;

const K = units.K;
const C = units.C;
const Pa = units.Pa;
const KPa = units.KPa;
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
const WATER_GAS_CONSTANT = constants.WATER_GAS_CONSTANT;
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

    pub fn init(liquid_frac: f64, vapor_frac: f64) SteamError!LiquidVapor {
        if (liquid_frac < 0.0 or vapor_frac < 0.0) return SteamError.InvalidPhaseFractions;
        const sum = liquid_frac + vapor_frac;
        if (!std.math.approxEqAbs(f64, sum, 1.0, 1e-6)) return SteamError.InvalidPhaseFractions;
        return LiquidVapor{
            .liquid_frac = liquid_frac,
            .vapor_frac = vapor_frac,
        };
    }
};

pub const CompositePhaseRegion = struct {
    liquid_vapor: LiquidVapor,

    pub fn liquidVapor(self: CompositePhaseRegion) LiquidVapor {
        return self.liquid_vapor;
    }
};

//pub const PtPoint = struct {
//pressure: Pressure,
//temperature: Temperature,
//};

//pub const SatQuery = union(enum) {
//SatTQuery: struct {
//temperature: Temperature,
//phase_region: SteamNonCriticalPhaseRegion,
//},
//SatPressureQuery: struct {
//pressure: Pressure,
//phase_region: SteamNonCriticalPhaseRegion,
//},
//};

//pub const SteamQuery = union(enum) {
//Pt: PtPoint,
//Sat: SatQuery,
//EntropyP: struct {
//entropy: EnergyPerMassTemperature,
//pressure: Pressure,
//},
//EnthalpyP: struct {
//enthalpy: EnergyPerMass,
//pressure: Pressure,
//},
//};

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
    pressure: Pa,
    temperature: K,
    tau: f64,
    pi: f64,
    gamma: f64,
    gamma_pi: f64,
    gamma_pi_pi: f64,
    gamma_tau: f64,
    gamma_tau_tau: f64,
    gamma_pi_tau: f64,
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

//fn extractPressure(query: SteamQuery) ?Pressure {
//return switch (query) {
//.Pt => |pt| pt.pressure,
//.Sat => |sat| switch (sat) {
//.SatPQuery => |sat_p| sat_p.pressure,
//else => null,
//},
//.EntropyP => |data| data.pressure,
//.EnthalpyP => |data| data.pressure,
//};
//}

//fn extractTemperature(query: SteamQuery) ?Temperature {
//return switch (query) {
//.Pt => |pt| pt.temperature,
//.Sat => |sat| switch (sat) {
//.SatTQuery => |sat_t| sat_t.temperature,
//else => null,
//},
//else => null,
//};
//}

fn checkIfOutOfRange(pressure: ?Pa, temperature: ?K) SteamError!void {
    const opt_p = if (pressure) |p| p.value else null;
    const opt_t = if (temperature) |t| t.value else null;

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

fn getRegionFromPressureAndTemperature(pressure: Pa, temperature: K) SteamError!Iapws97Region {
    const p = pressure.value;
    const t = temperature.value;

    if (t > 273.15 + 800.0) return Iapws97Region.Region5;
    if (t > 273.15 + 600.0) return Iapws97Region.Region2;

    const sat_pressure_optional = getSatPressure(temperature) catch |err| switch (err) {
        else => null,
    };

    if (sat_pressure_optional) |sat_pressure| {
        if (p == sat_pressure.value) return Iapws97Region.Region4;
        if (p < sat_pressure.value) return Iapws97Region.Region2;
        return Iapws97Region.Region1;
    }

    const boundary_pressure = getBoundary34Pressure(temperature) catch |err| switch (err) {
        else => return err,
    };

    if (p < boundary_pressure.value) return Iapws97Region.Region2;
    return Iapws97Region.Region3;
}

const SatRegion = struct { pressure: Pa, temperature: K, region: Iapws97Region };

fn getRegionFromSatTemperature(phase_region: SteamNonCriticalPhaseRegion, temperature: K) SteamError!SatRegion {
    return .{
        .pressure = try getSatPressure(temperature),
        .temperature = temperature,
        .region = if (phase_region == .Liquid) Iapws97Region.Region1 else Iapws97Region.Region2,
    };
}

fn getRegionFromSatPressure(phase_region: SteamNonCriticalPhaseRegion, pressure: Pa) SteamError!SatRegion {
    return .{
        .pressure = pressure,
        .temperature = try getSatTemperature(pressure),
        .region = if (phase_region == .Liquid) Iapws97Region.Region1 else Iapws97Region.Region2,
    };
}

fn createEntryFromRegionPoint(
    specific_region_point: SpecificRegionPoint,
    phase_region: PhaseKind,
) PtvEntry {
    const temperature = specific_region_point.temperature.value;
    const pressure = specific_region_point.pressure.value;
    const pi = specific_region_point.pi;
    const tau = specific_region_point.tau;
    const gamma = specific_region_point.gamma;
    const gamma_pi = specific_region_point.gamma_pi;
    const gamma_pi_pi = specific_region_point.gamma_pi_pi;
    const gamma_tau = specific_region_point.gamma_tau;
    const gamma_tau_tau = specific_region_point.gamma_tau_tau;
    const gamma_pi_tau = specific_region_point.gamma_pi_tau;

    const internal_energy = WATER_GAS_CONSTANT.value * temperature * (tau * gamma_tau - pi * gamma_pi);
    const enthalpy = WATER_GAS_CONSTANT.value * temperature * tau * gamma_tau;
    const entropy = WATER_GAS_CONSTANT.value * (tau * gamma_tau - gamma);
    const cv = WATER_GAS_CONSTANT.value *
        (-std.math.pow(f64, tau, @as(f64, 2)) * gamma_tau_tau +
            std.math.pow(f64, gamma_pi - tau * gamma_pi_tau, @as(f64, 2)) / gamma_pi_pi);
    const cp = WATER_GAS_CONSTANT.value * -std.math.pow(f64, tau, @as(f64, 2)) * gamma_tau_tau;
    const speed_of_sound = std.math.sqrt(WATER_GAS_CONSTANT.value * temperature *
        (std.math.pow(f64, gamma_pi, @as(f64, 2)) /
            ((std.math.pow(f64, gamma_pi - tau * gamma_pi_tau, @as(f64, 2)) /
                (std.math.pow(f64, tau, @as(f64, 2)) * gamma_tau_tau)) - gamma_pi_pi)));
    const specific_volume = pi * (gamma_pi * WATER_GAS_CONSTANT.value * temperature) / pressure;

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

fn gibbsMethod(pressure: Pa, temperature: K) PtvEntry {
    const pressure_value = pressure.value;
    const temperature_value = temperature.value;
    const pi = pressure_value / 16.53e6;
    const tau = 1386.0 / temperature_value;

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
        .pressure = pressure,
        .temperature = temperature,
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
    pressure: Pa,
    temperature: K,
    ideal_points: []const JnRegionPoint,
    residual_points: []const IjnRegionPoint,
) PtvEntry {
    const pressure_value = pressure.value;
    const temperature_value = temperature.value;
    const pi = pressure_value / 1.0e6;

    var gamma = @log(pi);
    var gamma_pi = @as(f64, 1) / pi;
    var gamma_pi_pi = @as(f64, -1) / std.math.pow(f64, pi, 2);
    var gamma_tau: f64 = 0;
    var gamma_tau_tau: f64 = 0;
    var gamma_pi_tau: f64 = 0;

    for (ideal_points) |region_point| {
        const n = region_point.n;
        const j = region_point.j;
        gamma += n * pow(f64, tau, j);
        gamma_tau += n * j * pow(f64, tau, j - @as(f64, 1));
        gamma_tau_tau += n * j * (j - @as(f64, 1)) * pow(f64, tau, j - @as(f64, 2));
    }

    for (residual_points) |region_point| {
        const n = region_point.n;
        const i = region_point.i;
        const j = region_point.j;
        const pow_pi = pow(f64, pi, i);
        const pow_tau = pow(f64, tau - tau_shift, j);
        gamma += n * pow_pi * pow_tau;
        gamma_pi += n * i * pow(f64, pi, i - @as(f64, 1)) * pow_tau;
        gamma_pi_pi += n * i * (i - @as(f64, 1)) * pow(f64, pi, i - @as(f64, 2)) * pow_tau;
        gamma_tau += n * pow_pi * j * pow(f64, tau - tau_shift, j - @as(f64, 1));
        gamma_tau_tau += n * pow_pi * j * (j - @as(f64, 1)) * pow(f64, tau - tau_shift, j - @as(f64, 2));
        gamma_pi_tau += n * i * pow(f64, pi, i - @as(f64, 1)) * j * pow(f64, tau - tau_shift, j - @as(f64, 1));
    }

    const phase_info = blk: {
        const above_critical_temp = temperature_value > constants.WATER_CRITICAL_TEMPERATURE.value;
        const above_critical_pressure = pressure_value > constants.WATER_CRITICAL_PRESSURE.value;
        if (above_critical_temp and above_critical_pressure) {
            break :blk PhaseKind{ .SupercriticalFluid = {} };
        } else if (above_critical_temp and !above_critical_pressure) {
            break :blk PhaseKind{ .Gas = {} };
        } else {
            break :blk PhaseKind{ .NonCritical = .Vapor };
        }
    };

    const specific_region_point = SpecificRegionPoint{
        .pressure = pressure,
        .temperature = temperature,
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

fn region3BySpecificVolume(temperature: K, specific_volume: f64) PtvEntry {
    const density = 1.0 / specific_volume;
    const n1 = steam_constants.REGION_3_N1.n;
    const delta = density / 322.0;
    const temperature_value = temperature.value;
    const tau = 647.096 / temperature_value;

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

    const pressure = phi_delta * delta * density * WATER_GAS_CONSTANT.value * temperature_value;
    const internal_energy = tau * phi_tau * WATER_GAS_CONSTANT.value * temperature_value;
    const enthalpy = (tau * phi_tau + delta * phi_delta) * WATER_GAS_CONSTANT.value * temperature_value;
    const entropy = (tau * phi_tau - phi) * WATER_GAS_CONSTANT.value;
    const cv = -std.math.pow(f64, tau, 2.0) * phi_tau_tau * WATER_GAS_CONSTANT.value;
    const cp = (-std.math.pow(f64, tau, 2.0) * phi_tau_tau +
        std.math.pow(f64, delta * phi_delta - delta * tau * phi_delta_tau, 2.0) /
            (2.0 * delta * phi_delta + std.math.pow(f64, delta, 2.0) * phi_delta_delta)) *
        WATER_GAS_CONSTANT.value;

    const speed_of_sound = std.math.sqrt(
        (2.0 * delta * phi_delta + std.math.pow(f64, delta, 2.0) * phi_delta_delta -
            std.math.pow(f64, delta * phi_delta - delta * tau * phi_delta_tau, 2.0) /
                (std.math.pow(f64, tau, 2.0) * phi_tau_tau)) * WATER_GAS_CONSTANT.value * temperature_value,
    );

    return PtvEntry{
        .temperature = Temperature{ .k = K.init(temperature_value) },
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

const Region3Tuple = struct { pressure: Pa, temperature: K };
fn region3Residual(region3_current_point: *Region3Tuple, x: f64) f64 {
    const entry = region3BySpecificVolume(region3_current_point.*.temperature, x);
    return entry.pressure.convertToSiUnit().value -
        region3_current_point.*.pressure.value;
}

fn region3Method(pressure: Pa, temperature: K) SteamError!PtvEntry {
    var c = root_finder.Closure(Region3Tuple){ .ctx = .{ .pressure = pressure, .temperature = temperature }, .func = region3Residual };
    const specific_volume = root_finder.secantMethod(Region3Tuple, &c, 1.0 / 500.0, 1e-4) catch {
        return SteamError.FailedToConverge;
    };
    return region3BySpecificVolume(temperature, specific_volume);
}

fn getEntryFromPressureAndTemperature(pressure: Pa, temperature: K, region: Iapws97Region) SteamError!PtvEntry {
    const temperature_value = temperature.value;
    return switch (region) {
        .Region1, .Region4 => gibbsMethod(pressure, temperature),
        .Region2 => vaporMethod(
            @as(f64, 540.0) / temperature_value,
            @as(f64, 0.5),
            pressure,
            temperature,
            steam_constants.REGION_2_IDEAL,
            steam_constants.REGION_2_RESIDUAL,
        ),
        .Region3 => try region3Method(pressure, temperature),
        .Region5 => vaporMethod(
            @as(f64, 1000) / temperature_value,
            @as(f64, 0),
            pressure,
            temperature,
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

    const phase_info = LiquidVapor.init(liq_frac, vap_frac) catch |err| {
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
    pressure_si: Pa,
    target_value: f64,
    comptime getPropValue: fn (PtvEntry) f64,
) SteamError!PtvEntry {
    const pressure = Pressure{ .pa = pressure_si };

    const liquid_entry_optional = getSteamEntryBySatPressure(.Liquid, pressure) catch null;
    const vapor_entry_optional = getSteamEntryBySatPressure(.Vapor, pressure) catch null;

    if (liquid_entry_optional) |liquid_entry| {
        if (vapor_entry_optional) |vapor_entry| {
            const liquid_value = getPropValue(liquid_entry);
            const vapor_value = getPropValue(vapor_entry);
            if (liquid_value <= target_value and vapor_value >= target_value) {
                const liq_frac = (vapor_value - target_value) /
                    (vapor_value - liquid_value);
                return interpolateEntry(liquid_entry, vapor_entry, liq_frac);
            }
        }
    }

    const ctx = struct {
        pressure: Pa,
        target: f64,

        fn eval(self: @This(), temperature: f64) f64 {
            const result = getSteamEntryByPressureAndTemperature(.{ .pa = self.pressure }, .{ .k = K.init(temperature) }) catch return std.math.nan(f64);
            return getPropValue(result) - self.target;
        }
    }{ .pressure = pressure_si, .target = target_value };

    const tol = 1e-5;
    if (tol < 0.0) return SteamError.FailedToConverge;
    const max_iter: usize = 50;
    const eps: f64 = 1e-4;

    var x0: f64 = 310.0;
    var y0 = ctx.eval(x0);
    if (!std.math.isFinite(y0)) return SteamError.FailedToConverge;
    if (@abs(y0) <= tol) {
        return getSteamEntryByPressureAndTemperature(.{ .pa = pressure_si }, .{ .k = K.init(x0) });
    }

    const x1_shift: f64 = x0 * (1.0 + eps);
    var x1: f64 = x1_shift + if (x1_shift >= 0.0) eps else -eps;
    var y1 = ctx.eval(x1);
    if (!std.math.isFinite(y1)) return SteamError.FailedToConverge;

    var tries: usize = 0;
    while (tries < max_iter) : (tries += 1) {
        if (@abs(y1) <= tol) {
            return getSteamEntryByPressureAndTemperature(.{ .pa = pressure_si }, .{ .k = K.init(x1) });
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
            return getSteamEntryByPressureAndTemperature(.{ .pa = pressure_si }, .{ .k = K.init(x2) });
        }

        x0 = x1;
        y0 = y1;
        x1 = x2;
        y1 = y2;
    }

    return SteamError.FailedToConverge;
}

pub fn getSteamEntryByPressureAndTemperature(pressure: Pressure, temperature: Temperature) SteamError!PtvEntry {
    const pressure_si = pressure.convertToSiUnit();
    const temperature_si = temperature.convertToSiUnit();

    try checkIfOutOfRange(pressure_si, temperature_si);

    const region = try getRegionFromPressureAndTemperature(pressure_si, temperature_si);

    return try getEntryFromPressureAndTemperature(pressure_si, temperature_si, region);
}

pub fn getSteamEntryBySatPressure(phase_region: SteamNonCriticalPhaseRegion, pressure: Pressure) SteamError!PtvEntry {
    const pressure_si = pressure.convertToSiUnit();

    try checkIfOutOfRange(pressure_si, null);

    const result = try getRegionFromSatPressure(phase_region, pressure_si);

    return try getEntryFromPressureAndTemperature(result.pressure, result.temperature, result.region);
}

pub fn getSteamEntryBySatTemperature(phase_region: SteamNonCriticalPhaseRegion, temperature: Temperature) SteamError!PtvEntry {
    const temperature_si = temperature.convertToSiUnit();

    try checkIfOutOfRange(null, temperature_si);

    const result = try getRegionFromSatTemperature(phase_region, temperature_si);

    return try getEntryFromPressureAndTemperature(result.pressure, result.temperature, result.region);
}

pub fn getSteamEntryByPressureAndEntropy(pressure: Pressure, entropy: EnergyPerMassTemperature) SteamError!PtvEntry {
    const pressure_si = pressure.convertToSiUnit();
    try checkIfOutOfRange(pressure_si, null);

    return try iteratePtEntrySolution(
        pressure_si,
        entropy.convertToSiUnit().value,
        struct {
            fn get(entry: PtvEntry) f64 {
                return entry.entropy.convertToSiUnit().value;
            }
        }.get,
    );
}

pub fn getSteamEntryByPressureAndEnthalpy(pressure: Pressure, enthalpy: EnergyPerMass) SteamError!PtvEntry {
    const pressure_si = pressure.convertToSiUnit();
    try checkIfOutOfRange(pressure_si, null);

    return try iteratePtEntrySolution(
        pressure_si,
        enthalpy.convertToSiUnit().value,
        struct {
            fn get(entry: PtvEntry) f64 {
                return entry.enthalpy.convertToSiUnit().value;
            }
        }.get,
    );
}

fn expectPtvPointsAreEqual(expected: PtvEntry, actual: PtvEntry) !void {
    std.debug.print("testing PtvEntry pressure\n", .{});
    try std.testing.expectApproxEqAbs(expected.pressure.convertToSiUnit().value, actual.pressure.convertToSiUnit().value, 1e-3);
    std.debug.print("testing PtvEntry temperature\n", .{});
    try std.testing.expectApproxEqAbs(expected.temperature.convertToSiUnit().value, actual.temperature.convertToSiUnit().value, 1e-3);
    std.debug.print("testing PtvEntry phase_region\n", .{});
    if (expected.phase_region == .Composite and actual.phase_region == .Composite) {
        std.debug.print("testing PtvEntry liquid composite\n", .{});
        try std.testing.expectApproxEqAbs(expected.phase_region.Composite.liquid_vapor.liquid_frac, expected.phase_region.Composite.liquid_vapor.liquid_frac, 1e-3);
        std.debug.print("testing PtvEntry vapor composite\n", .{});
        try std.testing.expectApproxEqAbs(expected.phase_region.Composite.liquid_vapor.vapor_frac, expected.phase_region.Composite.liquid_vapor.vapor_frac, 1e-3);
    } else {
        try std.testing.expectEqual(expected.phase_region, actual.phase_region);
    }
    std.debug.print("testing PtvEntry internal_energy\n", .{});
    try std.testing.expectApproxEqAbs(expected.internal_energy.convertToSiUnit().value, actual.internal_energy.convertToSiUnit().value, 1e-2);
    try std.testing.expectApproxEqAbs(expected.internal_energy.convertToSiUnit().value, actual.internal_energy.convertToSiUnit().value, 1e-2);
    std.debug.print("testing PtvEntry enthalpy\n", .{});
    try std.testing.expectApproxEqAbs(expected.enthalpy.convertToSiUnit().value, actual.enthalpy.convertToSiUnit().value, 1e-2);
    try std.testing.expectApproxEqAbs(expected.enthalpy.convertToSiUnit().value, actual.enthalpy.convertToSiUnit().value, 1e-2);
    std.debug.print("testing PtvEntry entropy\n", .{});
    try std.testing.expectApproxEqAbs(expected.entropy.convertToSiUnit().value, actual.entropy.convertToSiUnit().value, 1e-3);
    std.debug.print("testing PtvEntry cv\n", .{});
    try std.testing.expectApproxEqAbs(expected.cv.convertToSiUnit().value, actual.cv.convertToSiUnit().value, 1e-3);
    std.debug.print("testing PtvEntry cp\n", .{});
    try std.testing.expectApproxEqAbs(expected.cp.convertToSiUnit().value, actual.cp.convertToSiUnit().value, 1e-3);
    std.debug.print("testing PtvEntry speed_of_sound\n", .{});
    try std.testing.expectApproxEqAbs(expected.speed_of_sound.convertToSiUnit().value, actual.speed_of_sound.convertToSiUnit().value, 1e-3);
    std.debug.print("testing PtvEntry specific_volume\n", .{});
    try std.testing.expectApproxEqAbs(expected.specific_volume.convertToSiUnit().value, actual.specific_volume.convertToSiUnit().value, 1e-3);
}

test "Pressure and Temperature Query Test Region 3" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(750) },
        .pressure = Pressure{ .kpa = KPa.init(78.309563916917e3) },
        .phase_region = PhaseKind{ .SupercriticalFluid = {} },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(2102.069317626429e3) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(2258.688445460262e3) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(4.469719056217e3) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.71701677121e3) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(6.341653594791e3) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(760.696040876798) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(@as(f64, 1) / @as(f64, 500)) },
    };

    const actual = try getSteamEntryByPressureAndTemperature(expected.pressure, expected.temperature);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Pressure and Temperature Query Test Region 1" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(473.15) },
        .pressure = Pressure{ .pa = Pa.init(40e6) },
        .phase_region = PhaseKind{ .NonCritical = .Liquid },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(825.228016170348e3) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(870.124259682489e3) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.275752861241e3) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(3.292858637199e3) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(4.315767590903e3) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(1457.418351596083) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(0.001122406088) },
    };

    const actual = try getSteamEntryByPressureAndTemperature(expected.pressure, expected.temperature);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Pressure and Temperature Query Test Region 5" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(2000.0) },
        .pressure = Pressure{ .pa = Pa.init(30e6) },
        .phase_region = PhaseKind{ .SupercriticalFluid = {} },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(5637.070382521894e3) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(6571.226038618478e3) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(8.536405231138e3) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.395894362358e3) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.885698818781e3) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(1067.369478777425) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(0.03113852187) },
    };

    const actual = try getSteamEntryByPressureAndTemperature(expected.pressure, expected.temperature);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Pressure and Temperature Query Test Region 2" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(823.15) },
        .pressure = Pressure{ .pa = Pa.init(14e6) },
        .phase_region = PhaseKind{ .Gas = {} },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(3114.302136294585e3) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(3460.987255128561e3) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(6.564768889364e3) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(1.892708832325e3) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.666558503968e3) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(666.050616844223) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(0.024763222774) },
    };

    const actual = try getSteamEntryByPressureAndTemperature(expected.pressure, expected.temperature);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Saturated Pressure Liquid Query Test Region 1" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(393.361545936488) },
        .pressure = Pressure{ .pa = Pa.init(0.2e6) },
        .phase_region = PhaseKind{ .NonCritical = .Liquid },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(504471.741847973) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(504683.84552926) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(1530.0982011075) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(3666.99397284121) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(4246.73524917536) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(1520.69128792808) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(0.00106051840643552) },
    };

    const actual = try getSteamEntryBySatPressure(.Liquid, expected.pressure);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Saturated Pressure Vapor Query Test Region 2" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(393.361545936488) },
        .pressure = Pressure{ .pa = Pa.init(0.2e6) },
        .phase_region = PhaseKind{ .NonCritical = .Vapor },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(2529094.32835793) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(2706241.34137425) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(7126.8563914686) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(1615.96336473298) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2175.22318865273) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(481.883535821489) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(0.885735065081644) },
    };

    const actual = try getSteamEntryBySatPressure(.Vapor, expected.pressure);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Saturated Temperature Liquid Query Test Region 1" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(393.361545936488) },
        .pressure = Pressure{ .pa = Pa.init(0.2e6) },
        .phase_region = PhaseKind{ .NonCritical = .Liquid },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(504471.741847973) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(504683.84552926) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(1530.0982011075) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(3666.99397284121) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(4246.73524917536) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(1520.69128792808) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(0.00106051840643552) },
    };

    const actual = try getSteamEntryBySatTemperature(.Liquid, expected.temperature);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Saturated Temperature Vapor Query Test Region 2" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(393.361545936488) },
        .pressure = Pressure{ .pa = Pa.init(0.2e6) },
        .phase_region = PhaseKind{ .NonCritical = .Vapor },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(2529094.32835793) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(2706241.34137425) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(7126.8563914686) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(1615.96336473298) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2175.22318865273) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(481.883535821489) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(0.885735065081644) },
    };

    const actual = try getSteamEntryBySatTemperature(.Vapor, expected.temperature);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Entropy and Pressure Query Test Region 3" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(750) },
        .pressure = Pressure{ .kpa = KPa.init(78.309563916917e3) },
        .phase_region = PhaseKind{ .SupercriticalFluid = {} },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(2102.069317626429e3) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(2258.688445460262e3) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(4.469719056217e3) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.71701677121e3) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(6.341653594791e3) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(760.696040876798) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(@as(f64, 1) / @as(f64, 500)) },
    };

    const actual = try getSteamEntryByPressureAndEntropy(expected.pressure, expected.entropy);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Entropy and Pressure Query Test Region 5" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(2000.0) },
        .pressure = Pressure{ .pa = Pa.init(30e6) },
        .phase_region = PhaseKind{ .SupercriticalFluid = {} },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(5637.070382521894e3) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(6571.226038618478e3) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(8.536405231138e3) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.395894362358e3) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.885698818781e3) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(1067.369478777425) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(0.03113852187) },
    };

    const actual = try getSteamEntryByPressureAndEntropy(expected.pressure, expected.entropy);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Entropy and Pressure Query Test Region 2" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(823.15) },
        .pressure = Pressure{ .pa = Pa.init(14e6) },
        .phase_region = PhaseKind{ .Gas = {} },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(3114.302136294585e3) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(3460.987255128561e3) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(6.564768889364e3) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(1.892708832325e3) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.666558503968e3) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(666.050616844223) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(0.024763222774) },
    };

    const actual = try getSteamEntryByPressureAndEntropy(expected.pressure, expected.entropy);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Entropy and Pressure Query Test Composite" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(318.957548207023) },
        .pressure = Pressure{ .pa = Pa.init(10e3) },
        .phase_region = PhaseKind{
            .Composite = CompositePhaseRegion{
                .liquid_vapor = try LiquidVapor.init(0.1950875529218672, 1 - 0.1950875529218672),
            },
        },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(1999135.82661328) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(2117222.94886314) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(6.6858e3) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(1966.28009225455) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2377.86300751001) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(655.005141924186) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(1.0 / 193.16103883) },
    };

    const actual = try getSteamEntryByPressureAndEntropy(expected.pressure, expected.entropy);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Enthalpy and Pressure Query Test Region 3" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(750) },
        .pressure = Pressure{ .kpa = KPa.init(78.309563916917e3) },
        .phase_region = PhaseKind{ .SupercriticalFluid = {} },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(2102.069317626429e3) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(2258.688445460262e3) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(4.469719056217e3) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.71701677121e3) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(6.341653594791e3) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(760.696040876798) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(@as(f64, 1) / @as(f64, 500)) },
    };

    const actual = try getSteamEntryByPressureAndEnthalpy(expected.pressure, expected.enthalpy);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Enthalpy and Pressure Query Test Region 1" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(473.15) },
        .pressure = Pressure{ .pa = Pa.init(40e6) },
        .phase_region = PhaseKind{ .NonCritical = .Liquid },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(825.228016170348e3) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(870.124259682489e3) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.275752861241e3) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(3.292858637199e3) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(4.315767590903e3) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(1457.418351596083) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(0.001122406088) },
    };

    const actual = try getSteamEntryByPressureAndEnthalpy(expected.pressure, expected.enthalpy);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Enthalpy and Pressure Query Test Region 5" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(2000.0) },
        .pressure = Pressure{ .pa = Pa.init(30e6) },
        .phase_region = PhaseKind{ .SupercriticalFluid = {} },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(5637.070382521894e3) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(6571.226038618478e3) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(8.536405231138e3) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.395894362358e3) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.885698818781e3) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(1067.369478777425) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(0.03113852187) },
    };

    const actual = try getSteamEntryByPressureAndEnthalpy(expected.pressure, expected.enthalpy);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Enthalpy and Pressure Query Test Region 2" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(823.15) },
        .pressure = Pressure{ .pa = Pa.init(14e6) },
        .phase_region = PhaseKind{ .Gas = {} },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(3114.302136294585e3) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(3460.987255128561e3) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(6.564768889364e3) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(1.892708832325e3) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2.666558503968e3) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(666.050616844223) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(0.024763222774) },
    };

    const actual = try getSteamEntryByPressureAndEnthalpy(expected.pressure, expected.enthalpy);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Enthalpy and Pressure Query Test Composite" {
    const expected = PtvEntry{
        .temperature = Temperature{ .k = K.init(318.957548207023) },
        .pressure = Pressure{ .pa = Pa.init(10e3) },
        .phase_region = PhaseKind{
            .Composite = CompositePhaseRegion{
                .liquid_vapor = try LiquidVapor.init(0.1950875529218672, 1 - 0.1950875529218672),
            },
        },
        .internal_energy = EnergyPerMass{ .j_per_kg = JPerKg.init(1999135.82661328) },
        .enthalpy = EnergyPerMass{ .j_per_kg = JPerKg.init(2117222.94886314) },
        .entropy = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(6.6858e3) },
        .cv = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(1966.28009225455) },
        .cp = EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(2377.86300751001) },
        .speed_of_sound = Velocity{ .m_per_sec = MPerSec.init(655.005141924186) },
        .specific_volume = SpecificVolume{ .m3_per_kg = M3PerKg.init(1.0 / 193.16103883) },
    };

    const actual = try getSteamEntryByPressureAndEnthalpy(expected.pressure, expected.enthalpy);
    try expectPtvPointsAreEqual(expected, actual);
}

test "Low Temperature Error Test 1" {
    const temperature = Temperature{ .k = K.init(273.0) };
    const pressure = Pressure{ .pa = Pa.init(40e6) };

    const actual = getSteamEntryByPressureAndTemperature(pressure, temperature);
    try std.testing.expectError(SteamError.TemperatureLow, actual);
}

test "Low Temperature Error Test 2" {
    const temperature = Temperature{ .k = K.init(273.0) };
    const pressure = Pressure{ .pa = Pa.init(60e6) };

    const actual = getSteamEntryByPressureAndTemperature(pressure, temperature);
    try std.testing.expectError(SteamError.TemperatureLow, actual);
}

test "High Temperature Error Test 1" {
    const temperature = Temperature{ .c = C.init(2001) };
    const pressure = Pressure{ .pa = Pa.init(40e6) };

    const actual = getSteamEntryByPressureAndTemperature(pressure, temperature);
    try std.testing.expectError(SteamError.TemperatureHigh, actual);
}

test "High Temperature Error Test 2" {
    const temperature = Temperature{ .c = C.init(801) };
    const pressure = Pressure{ .pa = Pa.init(60e6) };

    const actual = getSteamEntryByPressureAndTemperature(pressure, temperature);
    try std.testing.expectError(SteamError.TemperatureHigh, actual);
}

test "Low Pressure Error Test 1" {
    const temperature = Temperature{ .c = C.init(799) };
    const pressure = Pressure{ .pa = Pa.init(-1) };

    const actual = getSteamEntryByPressureAndTemperature(pressure, temperature);
    try std.testing.expectError(SteamError.PressureLow, actual);
}

test "Low Pressure Error Test 2" {
    const temperature = Temperature{ .c = C.init(801) };
    const pressure = Pressure{ .pa = Pa.init(-1) };

    const actual = getSteamEntryByPressureAndTemperature(pressure, temperature);
    try std.testing.expectError(SteamError.PressureLow, actual);
}

test "High Temperature Error Test 3" {
    const temperature = Temperature{ .c = C.init(801) };
    const pressure = Pressure{ .pa = Pa.init(51e6) };

    const actual = getSteamEntryByPressureAndTemperature(pressure, temperature);
    try std.testing.expectError(SteamError.TemperatureHigh, actual);
}

test "High Pressure Error Test 3" {
    const temperature = Temperature{ .c = C.init(799) };
    const pressure = Pressure{ .pa = Pa.init(101e6) };

    const actual = getSteamEntryByPressureAndTemperature(pressure, temperature);
    try std.testing.expectError(SteamError.PressureHigh, actual);
}
