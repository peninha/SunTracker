Date.prototype.isLeapYear = function() {
    var year = this.getFullYear();
    if((year & 3) != 0) return false;
    return ((year % 100) != 0 || (year % 400) == 0);
};

// Get Day of Year
Date.prototype.getDOY = function() {
    var dayCount = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    var mn = this.getUTCMonth();
    var dn = this.getUTCDate();
    var dayOfYear = dayCount[mn] + dn;
    if(mn > 1 && this.isLeapYear()) dayOfYear++;
    return dayOfYear;
};

(function() {
RIGHT_ASCENSION = 14.2612 // hours
// DECLINATION = 19.1873 // degrees
DECLINATION = 0.33488 // radians

// get latitude and longitude from ... something 
const LAT = 37.338;
const LNG = -121.886;

GST_2016 = 6.821 // base for jan 01 0h UT

const now = new Date();
const dayOfYear = now.getDOY();
const ut_in_hours = now.getUTCHours() + now.getUTCMinutes()/60 + now.getUTCSeconds()/(60*60) + now.getUTCMilliseconds()/(60*60*1000);
cur_gst = (GST_2016 + 0.06570982441908*dayOfYear + 1.00273790935*ut_in_hours) % 24;

// 15 deg is 1 hr
// ensure above 0?
local_sidereal_time = (cur_gst + (LNG / 15) + 24) % 24;

console.log(`Local sidereal time (hours) ${local_sidereal_time}`);
console.log(`Local sidereal time (radians) ${toRadians(local_sidereal_time * 15)}`);

local_hour_angle = (local_sidereal_time - RIGHT_ASCENSION + 24) % 24
const localHourRadian = toRadians(local_hour_angle * 15);

function toRadians(degrees) {
	return degrees * Math.PI / 180;
}

function toDegrees(radians) {
    return radians * 180 / Math.PI;
}

const latRad = toRadians(LAT);
const decl = DECLINATION;
const adjacent = Math.cos(localHourRadian)*Math.sin(latRad) - Math.tan(decl) * Math.cos(latRad);
azimuth = Math.atan2(Math.sin(localHourRadian), adjacent);
console.log(`azimuth (radians) ${azimuth}`);
console.log(`azimuth (degrees) ${toDegrees(azimuth)}`);

altitude = Math.asin(Math.sin(latRad)*Math.sin(decl) + Math.cos(latRad)*Math.cos(decl)*Math.cos(localHourRadian));
console.log(`altitude (radians) ${altitude}`);
console.log(`altitude (degrees) ${toDegrees(altitude)}`);
})();