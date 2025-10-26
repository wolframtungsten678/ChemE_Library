const common = @import("common.zig");

const UnitValueType = common.UnitValueType;
const UnitLabel = common.UnitLabel;
const RawUnit = common.RawUnit;
const ParseUnitError = common.ParseUnitError;
const isMatch = common.isMatch;
const pascals_per_kpa = common.pascals_per_kpa;
const pascals_per_psi = common.pascals_per_psi;

pub const Pa = UnitValueType();
pub const KPa = UnitValueType();
pub const Lbf = UnitValueType();

const pressure_labels = [_]UnitLabel{
    .{ .abbreviation = "Pa", .plural = "pascals" },
    .{ .abbreviation = "kPa", .plural = "kilopascals" },
    .{ .abbreviation = "lbf/in²", .plural = "pounds-force per square inch" },
};

pub const Pressure = union(enum) {
    pa: Pa,
    kpa: KPa,
    lbf: Lbf,

    pub fn convertToSiUnit(self: Pressure) Pa {
        return switch (self) {
            .pa => |val| val,
            .kpa => |val| Pa.init(val.value * pascals_per_kpa),
            .lbf => |val| Pa.init(val.value * pascals_per_psi),
        };
    }

    pub fn listUnitLabels() []const UnitLabel {
        return &pressure_labels;
    }

    pub fn getSiUnitLabel() UnitLabel {
        return pressure_labels[0];
    }

    pub fn getValue(self: Pressure) f64 {
        return switch (self) {
            .pa => |val| val.value,
            .kpa => |val| val.value,
            .lbf => |val| val.value,
        };
    }

    pub fn tryConvert(self: Pressure, unit_display: []const u8) ParseUnitError!Pressure {
        const value_si = self.convertToSiUnit().value;

        if (isMatch(unit_display, "Pa")) return .{ .pa = Pa.init(value_si) };
        if (isMatch(unit_display, "kPa")) return .{ .kpa = KPa.init(value_si / pascals_per_kpa) };
        if (isMatch(unit_display, "lbf/in²")) return .{ .lbf = Lbf.init(value_si / pascals_per_psi) };

        return ParseUnitError.UnknownUnit;
    }

    pub fn fromRawUnit(raw: RawUnit) ParseUnitError!Pressure {
        if (isMatch(raw.unit_display, "Pa")) return .{ .pa = Pa.init(raw.value) };
        if (isMatch(raw.unit_display, "kPa")) return .{ .kpa = KPa.init(raw.value) };
        if (isMatch(raw.unit_display, "lbf/in²")) return .{ .lbf = Lbf.init(raw.value) };
        return ParseUnitError.UnknownUnit;
    }

    pub fn toRawUnit(self: Pressure) RawUnit {
        return switch (self) {
            .pa => |val| .{ .value = val.value, .unit_display = "Pa" },
            .kpa => |val| .{ .value = val.value, .unit_display = "kPa" },
            .lbf => |val| .{ .value = val.value, .unit_display = "lbf/in²" },
        };
    }
};
