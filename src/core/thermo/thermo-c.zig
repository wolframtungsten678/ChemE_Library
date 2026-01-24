const iapws97 = @import("steam/iapws97.zig");
const units = @import("../units.zig");
const Pressure = units.Pressure;
const Pa = units.Pa;
const Temperature = units.Temperature;
const K = units.K;

pub const CSteamResult = extern struct {
    ok: bool,
    err_code: c_int,
    phase_kind: c_int,
    phase_region: c_int,
    liquid_frac: f64,
    vapor_frac: f64,
    pressure: f64,
    temperature: f64,
    internal_energy: f64,
    enthalpy: f64,
    entropy: f64,
    cv: f64,
    cp: f64,
    speed_of_sound: f64,
    specific_volume: f64,
};

const PhaseKindTag = struct {
    const SupercriticalFluid: c_int = 0;
    const Gas: c_int = 1;
    const NonCritical: c_int = 2;
    const Composite: c_int = 3;
    const Unknown: c_int = -1;
};

const NonCriticalTag = struct {
    const Liquid: c_int = 0;
    const Vapor: c_int = 1;
    const Unknown: c_int = -1;
};

fn convertToCSteamResult(res: iapws97.SteamError!iapws97.PtvEntry) CSteamResult {
    if (res) |ok| {
        var phase_kind: c_int = PhaseKindTag.Unknown;
        var phase_region: c_int = NonCriticalTag.Unknown;
        var liquid_frac: f64 = 0;
        var vapor_frac: f64 = 0;

        switch (ok.phase_region) {
            .SupercriticalFluid => phase_kind = PhaseKindTag.SupercriticalFluid,
            .Gas => phase_kind = PhaseKindTag.Gas,
            .NonCritical => |region| {
                phase_kind = PhaseKindTag.NonCritical;
                phase_region = switch (region) {
                    .Liquid => NonCriticalTag.Liquid,
                    .Vapor => NonCriticalTag.Vapor,
                };
            },
            .Composite => |composite| {
                phase_kind = PhaseKindTag.Composite;
                liquid_frac = composite.liquid_vapor.liquid_frac;
                vapor_frac = composite.liquid_vapor.vapor_frac;
            },
        }

        return .{
            .ok = true,
            .err_code = 0,
            .phase_kind = phase_kind,
            .phase_region = phase_region,
            .liquid_frac = liquid_frac,
            .vapor_frac = vapor_frac,
            .pressure = ok.pressure.getValue(),
            .temperature = ok.temperature.getValue(),
            .internal_energy = ok.internal_energy.getValue(),
            .enthalpy = ok.enthalpy.getValue(),
            .entropy = ok.entropy.getValue(),
            .cv = ok.cv.getValue(),
            .cp = ok.cp.getValue(),
            .speed_of_sound = ok.speed_of_sound.getValue(),
            .specific_volume = ok.specific_volume.getValue(),
        };
    } else |err| {
        return .{
            .ok = false,
            .err_code = @intFromError(err),
            .phase_kind = PhaseKindTag.Unknown,
            .phase_region = NonCriticalTag.Unknown,
            .liquid_frac = 0,
            .vapor_frac = 0,
            .pressure = 0,
            .temperature = 0,
            .internal_energy = 0,
            .enthalpy = 0,
            .entropy = 0,
            .cv = 0,
            .cp = 0,
            .speed_of_sound = 0,
            .specific_volume = 0,
        };
    }
}

pub export fn getSteamEntryByPressureAndTemperature(raw_p: f64, raw_t: f64) CSteamResult {
    const pressure = Pressure{ .pa = Pa.init(raw_p) };
    const temperature = Temperature{ .k = K.init(raw_t) };
    const result = iapws97.getSteamEntryByPressureAndTemperature(pressure, temperature);
    return convertToCSteamResult(result);
}
