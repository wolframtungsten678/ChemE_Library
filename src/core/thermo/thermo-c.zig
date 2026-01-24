const std = @import("std");
const iapws97 = @import("steam/iapws97.zig");
const units = @import("../units.zig");
const Pressure = units.Pressure;
const Pa = units.Pa;
const Temperature = units.Temperature;
const K = units.K;

pub const CSteamResult = extern struct {
    ok: bool,
    err_code: c_int,
    pressure: f64,
    temperature: f64,
};

fn convertToCSteamResult(res: iapws97.SteamError!iapws97.PtvEntry) CSteamResult {
    return if (res) |ok| .{
        .ok = true,
        .err_code = 0,
        .pressure = ok.pressure.getValue(),
        .temperature = ok.temperature.getValue(),
    } else |err| .{
        .ok = false,
        .err_code = @intFromError(err),
        .pressure = 0,
        .temperature = 0,
    };
}

pub export fn getSteamEntryByPressureAndTemperature(raw_p: f64, raw_t: f64) CSteamResult {
    const pressure = Pressure{ .pa = Pa.init(raw_p) };
    const temperature = Temperature{ .k = K.init(raw_t) };
    const result = iapws97.getSteamEntryByPressureAndTemperature(pressure, temperature);
    return convertToCSteamResult(result);
}
