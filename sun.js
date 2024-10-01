// convert Gregorian date to Julian day
// Ref: Astronomical Algorithms by Jean Meeus
function julianDay(year, month, day, hour, min, sec){ //Valid for Gregorian dates after 1582-10-15
	hour = hour || 0;
	min = min || 0;
	sec = sec || 0;
	if (month <= 2) {
		year -= 1;
		month += 12;
	};
	day += (hour)/24.0 + min/1440.0 + sec/86400.0;
	let A = Math.floor(year/100);
	let B = 2 - A + Math.floor(A/4);
	if (year < 1582 || (year==1582 && (month < 10 || (month==10 && day < 15)))){
		B = 0;
	}
	let JD = Math.floor(365.25*(year+4716)) + Math.floor(30.6001*(month + 1)) + day + B - 1524.5;
	return JD;
}

function JD2Date(JD,format){
	format = format || "day";
	JD += 0.5
	let Z = Math.floor(JD);
	let F = JD - Z;
	let A = Z;
	if (Z >= 2299161) {
		var alpha = Math.floor((Z - 1867216.25)/36524.25);
		A = Z + 1 + alpha - Math.floor(alpha/4);
	}
	let B = A + 1524;
	let C = Math.floor((B - 122.1)/365.25);
	let D = Math.floor(365.25 * C);
	let E = Math.floor((B - D)/30.6001);
	let day = B - D - Math.floor(30.6001 * E) + F;
	let month = E < 14 ? E - 1 : E - 13;
	let year = month > 2 ? C - 4716 : C - 4715;
	if (format == "day"){
		return {year : year, month : month, day : day};
	}else{
		let time = dayDec2Hour(day);
		let date = new Date(year, month-1, Math.floor(day), time.hour, time.min, Math.floor(time.sec), (time.sec - Math.floor(time.sec))*1000);
		let Y = new Date(year,0,1);
		
		let decimal = (date.getTime() - Y.getTime())/(24*60*60*numberDaysInYear(year)*1000);
		return year + decimal;
	}	
}

function numberDaysInYear(year){
	return isLeapYear(year) ? 366 : 365;
}

function isLeapYear(year){
	return ((year % 4 == 0) && (year % 100 != 0)) || (year % 400 == 0);
}

function dayDec2Hour(day){
	day = new Decimal(day);
	let hour = day.minus(day.floor()).times(24);
	let min = hour.minus(hour.floor()).times(60);
	let sec = min.minus(min.floor()).times(60);
	return {
		hour : hour.floor(),
		min : min.floor(),
		sec : sec.toFixed(4)
	};
}

// convert Equatorial System to Horizontal System
// Ref: https://en.wikipedia.org/wiki/Celestial_coordinate_system#Equatorial_%E2%86%94_horizontal
function sunHorizontPosition(jd, coords, UTC){
	var D = jd - 2451545.0;
	var g = DMath.fixAngle(357.529 + 0.98560028* D);
	var q = DMath.fixAngle(280.459 + 0.98564736* D);
	var L = DMath.fixAngle(q + 1.915* DMath.sin(g) + 0.020* DMath.sin(2*g));

	var R = 1.00014 - 0.01671*DMath.cos(g) - 0.00014*DMath.cos(2*g);
	var e = 23.439 - 0.00000036* D;

	var RA = DMath.arctan2(DMath.cos(e)* DMath.sin(L), DMath.cos(L))/ 15;
	var eqt = q/15 - DMath.fixHour(RA);
	var Decl = DMath.arcsin(DMath.sin(e)* DMath.sin(L));
	
	var jd0 = Math.trunc(jd - 0.5) + 0.5;
	var D0 = jd0 - 2451545.0;
	var H = (jd - jd0)*24;
	var T = D/36525.0;
	var GMST = (6.697374558 + 0.06570982441908*D0 + 1.00273790935*H + 0.000026*T*T) % 24;
	var Omega = 125.04 - 0.052954*D;
	var L = 280.47 + 0.98565*D
	var Epsilon = 23.4393 - 0.0000004*D;
	var eqeq = (-0.000319*DMath.sin(Omega) - 0.000024*DMath.sin(2*L)) * DMath.cos(Epsilon);
	var GAST = GMST + eqeq;
	lng = 1* coords[0];
	lat = 1* coords[1];
	elv = coords[2] ? 1* coords[2] : 0;
	var LST = GAST + lng/15.0;
	var HourAngle = DMath.fixAngle((LST - RA)*15);
	
	var Azimuth = -DMath.arctan2(DMath.cos(Decl)*DMath.sin(HourAngle), -DMath.sin(lat)*DMath.cos(Decl)*DMath.cos(HourAngle) + DMath.cos(lat)*DMath.sin(Decl));
	Azimuth = Azimuth < 0 ? Azimuth + 360 : Azimuth;
	var Altitude = DMath.arcsin(DMath.sin(lat)*DMath.sin(Decl) + DMath.cos(lat)*DMath.cos(Decl)*DMath.cos(HourAngle));
	
	var Refraction = 1.02/DMath.tan(Altitude + 10.3/(Altitude + 5.11));
	//Altitude += Refraction/60;
	
	return{
		RA: RA,
		Declination: Decl,
		Azimuth: Azimuth,
		Altitude: Altitude,
		HourAngle: HourAngle,
		LST: LST,
	};
}

function sunPosition(JD, coords){
	let JDE = UT2TT(JD);
	//JDE = 2448908.5;
	

	//L = DMath.fixAngle(L*180/Math.PI + 180); // Convert to Degree and sum 180°
	//B = DMath.fixAngle(B*180/Math.PI*-1); // Convert to Degree and invert signal
	//R = R*1;
	
	//Convert to the FK5 System
	//let time = (JDE - 2451545)/365250;
	//let T = time*10;
	let T = (JDE - 2451545)/36525;
	//let Lambda = L.plus(T.times(-1.397)).plus(T.pow(2).times(-0.00031)).times(Math.PI).div(180); // L - 1°.397T - 0°.00031T² converted to Radians
	//let Lambda = L - 1.397*T - 0.00031*T*T; // L - 1°.397T - 0°.00031T²
	//L = L - 0.09033/3600; // -0".09033
	//B = B + 0.03916/3600*(DMath.cos(Lambda) - DMath.sin(Lambda)); // +0".03916(cos(Lamda) - sin(Lambda))
	
	//Nutation correction
	 // Julian centuries from Epoch J2000.0
	//let Omega = T.pow(3).div(450000).plus(T.pow(2).times(0.0020708)).plus(T.times(-1934.136261)).plus(125.04452); // Longitude of the Ascending node of the Moon's mean orbit on the ecliptic, measured from the mean equinox of the date = 125.04452 - 1934.136261T + 0.0020708T² + T³/450000
	//T = time * 10;
	
	//let Omega = T.times(-1934.136261).plus(125.04452); // Longitude of the Ascending node of the Moon's mean orbit on the ecliptic, measured from the mean equinox of the date = 125.04452 - 1934.136261T + 0.0020708T² + T³/450000
	let Omega = DMath.fixAngle(125.04452 - 1934.136261*T + 0.0020708*T*T + T*T*T/450000); // Longitude of the Ascending node of the Moon's mean orbit on the ecliptic, measured from the mean equinox of the date = 125.04452 - 1934.136261T + 0.0020708T² + T³/450000
	
	//let Lmean = T.times(36000.7698).plus(280.4665); // Mean Longitude of the Sun = 280°.4665 + 36000°.7698T
	let Lmean = 280.4665 + 36000.7698*T; // Mean Longitude of the Sun = 280°.4665 + 36000°.7698T
	
	//let Lprimemean = T.times(481267.8813).plus(218.3165); // Mean Longitude of the Moon = 218°.3165 + 481267°.8813T
	let Lprimemean = 218.3165 + 481267.8813*T; // Mean Longitude of the Moon = 218°.3165 + 481267°.8813T
	
	//let DeltaPsi = Omega.sin().times(-17.20/3600).plus(Lmean.times(2).sin().times(-1.32/3600)).plus(Lprimemean.times(2).sin().times(-0.23/3600)).plus(Omega.times(2).sin().times(0.21/3600)); // Nutation in longitude = -17".20sin(Omega) -1".32sin(2Lmean) -0".23sin(2Lprimemean) +0".21sin(2Omega)
	let DeltaPsi = -17.20/3600*DMath.sin(Omega) - 1.32/3600*DMath.sin(2*Lmean) - 0.23/3600*DMath.sin(2*Lprimemean) + 0.21/3600*DMath.sin(2*Omega); // Nutation in longitude = -17".20sin(Omega) -1".32sin(2Lmean) -0".23sin(2Lprimemean) +0".21sin(2Omega)
	
	//let DeltaEpsilon = Omega.cos().times(9.20/3600).plus(Lmean.times(2).cos().times(0.57/3600)).plus(Lprimemean.times(2).cos().times(0.10/3600)).plus(Omega.times(2).cos().times(-0.09/3600)); // Nutation in obliquity = 9".20cos(Omega) +0".57cos(2Lmean) +0".10cos(2Lprimemean) -0".09cos(2Omega)
	let DeltaEpsilon = 9.20/3600*DMath.cos(Omega) + 0.57/3600*DMath.cos(2*Lmean) + 0.10/3600*DMath.cos(2*Lprimemean) - 0.09/3600*DMath.cos(2*Omega); // Nutation in obliquity = 9".20cos(Omega) +0".57cos(2Lmean) +0".10cos(2Lprimemean) -0".09cos(2Omega)
	
	//let EpsilonZero = T.pow(3).times(0.001813/3600).plus(T.pow(2).times(-0.00059/3600)).plus(T.times(-46.8150/3600)).plus(23 + 26/60 + 21.448/3600); // Mean obliquity = 23°26'21".448 - 46".8150T - 0".00059T² + 0".001813T³
	let EpsilonZero = (23+26/60+21.448/3600) - 46.8150/3600*T - 0.00059/3600*T*T + 0.001813/3600*T*T*T; // Mean obliquity = 23°26'21".448 - 46".8150T - 0".00059T² + 0".001813T³
	
	//let Epsilon = EpsilonZero.plus(DeltaEpsilon);
	let Epsilon = EpsilonZero + DeltaEpsilon;
	//L += DeltaPsi;
	
	//Aberration correction
	//L += -20.4898/(3600*R); // -20".4898/R
	
	//let Alpha = DMath.arctan2(DMath.sin(L)*DMath.cos(Epsilon) - DMath.arctan(B)*DMath.sin(Epsilon), DMath.cos(L));
	//let Delta = DMath.arcsin(DMath.sin(B)*DMath.cos(Epsilon) + DMath.cos(B)*DMath.sin(Epsilon)*DMath.sin(L));
	
	JD = JDE;
	//let time = new Decimal((JDE - 2451545)/365250.0);
	
	var D = JDE - 2451545.0;
	var g = DMath.fixAngle(357.529 + 0.98560028*D);
	var q = DMath.fixAngle(280.459 + 0.98564736*D);
	var L = DMath.fixAngle(q + 1.915*DMath.sin(g) + 0.020*DMath.sin(2*g));

	var R = 1.00014 - 0.01671* DMath.cos(g) - 0.00014* DMath.cos(2*g);
	var e = 23.439 - 0.00000036* D;

	var RA = DMath.arctan2(DMath.cos(e)*DMath.sin(L), DMath.cos(L))/ 15;
	var decl = DMath.arcsin(DMath.sin(e)*DMath.sin(L));

	var RA2 = DMath.arctan2(DMath.cos(Epsilon)*DMath.sin(L), DMath.cos(L))/ 15;
	var decl2 = DMath.arcsin(DMath.sin(Epsilon)*DMath.sin(L));
	
	var RA3 = DMath.arctan2(DMath.sin(L)*DMath.cos(Epsilon), DMath.cos(L));
	var decl3 = DMath.arcsin(DMath.sin(Epsilon)*DMath.sin(L));
	
	var jd0 = Math.trunc(JD - 0.5) + 0.5;
	var D0 = jd0 - 2451545.0;
	var H = (JD - jd0)*24;
	T = D/36525.0;
	var GMST = (6.697374558 + 0.06570982441908*D0 + 1.00273790935*H + 0.000026*T*T) % 24;
	Omega = 125.04 - 0.052954*D;
	var L = 280.47 + 0.98565*D
	Epsilon = 23.4393 - 0.0000004*D;
	var eqeq = (-0.000319*DMath.sin(Omega) - 0.000024*DMath.sin(2*L)) * DMath.cos(Epsilon);
	var GAST = GMST + eqeq;
	lat = 1* coords[0];
	lng = 1* coords[1];
	elv = coords[2] ? 1* coords[2] : 0;
	var LST = GAST + lng/15.0;
	var h = (LST - RA)*15;
	
	var A = -DMath.arctan2(DMath.cos(decl)*DMath.sin(h), -DMath.sin(lat)*DMath.cos(decl)*DMath.cos(h) + DMath.cos(lat)*DMath.sin(decl));
	A = A < 0 ? A + 360 : A;
	var a = DMath.arcsin(DMath.sin(lat)*DMath.sin(decl) + DMath.cos(lat)*DMath.cos(decl)*DMath.cos(h));
	
	//return {declination: decl, equation: eqt};
	
	console.log(L+0.0);
	//L
	//B
	//R
	// VSOP87
	// RA = 198.37812
	// Decl = 7.78381
}

function degreeDec2Arc(angle){
	degree = Math.trunc(angle);
	min = Math.trunc((angle - degree)*60);
	sec = ((angle - degree)*60 - min)*60;
	return { degree: degree, min: min, sec: sec };
}
function degreeDec2Hour(angle){
	hour = Math.trunc(angle/15);
	min = Math.trunc((angle/15 - hour)*60);
	sec = ((angle/15 - hour)*60 - min)*60;
	return { hour: hour, min: min, sec: sec };
}

//---------------------- Degree-Based Math Class -----------------------
var DMath = {

	dtr: function(d) { return (d * Math.PI) / 180.0; },
	rtd: function(r) { return (r * 180.0) / Math.PI; },

	sin: function(d) { return Math.sin(this.dtr(d)); },
	cos: function(d) { return Math.cos(this.dtr(d)); },
	tan: function(d) { return Math.tan(this.dtr(d)); },

	arcsin: function(d) { return this.rtd(Math.asin(d)); },
	arccos: function(d) { return this.rtd(Math.acos(d)); },
	arctan: function(d) { return this.rtd(Math.atan(d)); },

	arccot: function(x) { return this.rtd(Math.atan(1/x)); },
	arctan2: function(y, x) { return this.rtd(Math.atan2(y, x)); },

	fixAngle: function(a) { return this.fix(a, 360); },
	fixHour:  function(a) { return this.fix(a, 24 ); },

	fix: function(a, b) {
		a = a- b* (Math.floor(a/ b));
		return (a < 0) ? a+ b : a;
	}
}

function UT2TT(JD) {
	year = JD2Date(JD,"year");
	
	//http://hemel.waarnemen.com/Computing/data/deltat.dat
	//http://maia.usno.navy.mil/ser7/deltat.data
	//http://maia.usno.navy.mil/ser7/deltat.preds
	let deltaTable = [
		{"year" : -700 , "deltaT" : 20400},
		{"year" : -600 , "deltaT" : 18800},
		{"year" : -500 , "deltaT" : 17190},
		{"year" : -400 , "deltaT" : 15530},
		{"year" : -300 , "deltaT" : 14080},
		{"year" : -200 , "deltaT" : 12790},
		{"year" : -100 , "deltaT" : 11640},
		{"year" : 0    , "deltaT" : 10580},
		{"year" : 100  , "deltaT" : 9600},
		{"year" : 200  , "deltaT" : 8640},
		{"year" : 300  , "deltaT" : 7680},
		{"year" : 400  , "deltaT" : 6700},
		{"year" : 500  , "deltaT" : 5710},
		{"year" : 600  , "deltaT" : 4740},
		{"year" : 700  , "deltaT" : 3810},
		{"year" : 800  , "deltaT" : 2960},
		{"year" : 900  , "deltaT" : 2200},
		{"year" : 1000 , "deltaT" : 1570},
		{"year" : 1100 , "deltaT" : 1090},
		{"year" : 1200 , "deltaT" : 740},
		{"year" : 1300 , "deltaT" : 490},
		{"year" : 1400 , "deltaT" : 320},
		{"year" : 1500 , "deltaT" : 200},
		{"year" : 1600 , "deltaT" : 120},
		{"year" : 1620 , "deltaT" : 124},
		{"year" : 1621 , "deltaT" : 119},
		{"year" : 1622 , "deltaT" : 115},
		{"year" : 1623 , "deltaT" : 110},
		{"year" : 1624 , "deltaT" : 106},
		{"year" : 1625 , "deltaT" : 102},
		{"year" : 1626 , "deltaT" : 98},
		{"year" : 1627 , "deltaT" : 95},
		{"year" : 1628 , "deltaT" : 91},
		{"year" : 1629 , "deltaT" : 88},
		{"year" : 1630 , "deltaT" : 85},
		{"year" : 1631 , "deltaT" : 82},
		{"year" : 1632 , "deltaT" : 79},
		{"year" : 1633 , "deltaT" : 77},
		{"year" : 1634 , "deltaT" : 74},
		{"year" : 1635 , "deltaT" : 72},
		{"year" : 1636 , "deltaT" : 70},
		{"year" : 1637 , "deltaT" : 67},
		{"year" : 1638 , "deltaT" : 65},
		{"year" : 1639 , "deltaT" : 63},
		{"year" : 1640 , "deltaT" : 62},
		{"year" : 1641 , "deltaT" : 60},
		{"year" : 1642 , "deltaT" : 58},
		{"year" : 1643 , "deltaT" : 57},
		{"year" : 1644 , "deltaT" : 55},
		{"year" : 1645 , "deltaT" : 54},
		{"year" : 1646 , "deltaT" : 53},
		{"year" : 1647 , "deltaT" : 51},
		{"year" : 1648 , "deltaT" : 50},
		{"year" : 1649 , "deltaT" : 49},
		{"year" : 1650 , "deltaT" : 48},
		{"year" : 1651 , "deltaT" : 47},
		{"year" : 1652 , "deltaT" : 46},
		{"year" : 1653 , "deltaT" : 45},
		{"year" : 1654 , "deltaT" : 44},
		{"year" : 1655 , "deltaT" : 43},
		{"year" : 1656 , "deltaT" : 42},
		{"year" : 1657 , "deltaT" : 41},
		{"year" : 1658 , "deltaT" : 40},
		{"year" : 1659 , "deltaT" : 38},
		{"year" : 1660 , "deltaT" : 37},
		{"year" : 1661 , "deltaT" : 36},
		{"year" : 1662 , "deltaT" : 35},
		{"year" : 1663 , "deltaT" : 34},
		{"year" : 1664 , "deltaT" : 33},
		{"year" : 1665 , "deltaT" : 32},
		{"year" : 1666 , "deltaT" : 31},
		{"year" : 1667 , "deltaT" : 30},
		{"year" : 1668 , "deltaT" : 28},
		{"year" : 1669 , "deltaT" : 27},
		{"year" : 1670 , "deltaT" : 26},
		{"year" : 1671 , "deltaT" : 25},
		{"year" : 1672 , "deltaT" : 24},
		{"year" : 1673 , "deltaT" : 23},
		{"year" : 1674 , "deltaT" : 22},
		{"year" : 1675 , "deltaT" : 21},
		{"year" : 1676 , "deltaT" : 20},
		{"year" : 1677 , "deltaT" : 19},
		{"year" : 1678 , "deltaT" : 18},
		{"year" : 1679 , "deltaT" : 17},
		{"year" : 1680 , "deltaT" : 16},
		{"year" : 1681 , "deltaT" : 15},
		{"year" : 1682 , "deltaT" : 14},
		{"year" : 1683 , "deltaT" : 14},
		{"year" : 1684 , "deltaT" : 13},
		{"year" : 1685 , "deltaT" : 12},
		{"year" : 1686 , "deltaT" : 12},
		{"year" : 1687 , "deltaT" : 11},
		{"year" : 1688 , "deltaT" : 11},
		{"year" : 1689 , "deltaT" : 10},
		{"year" : 1690 , "deltaT" : 10},
		{"year" : 1691 , "deltaT" : 10},
		{"year" : 1692 , "deltaT" : 9},
		{"year" : 1693 , "deltaT" : 9},
		{"year" : 1694 , "deltaT" : 9},
		{"year" : 1695 , "deltaT" : 9},
		{"year" : 1696 , "deltaT" : 9},
		{"year" : 1697 , "deltaT" : 9},
		{"year" : 1698 , "deltaT" : 9},
		{"year" : 1699 , "deltaT" : 9},
		{"year" : 1700 , "deltaT" : 9},
		{"year" : 1701 , "deltaT" : 9},
		{"year" : 1702 , "deltaT" : 9},
		{"year" : 1703 , "deltaT" : 9},
		{"year" : 1704 , "deltaT" : 9},
		{"year" : 1705 , "deltaT" : 9},
		{"year" : 1706 , "deltaT" : 9},
		{"year" : 1707 , "deltaT" : 9},
		{"year" : 1708 , "deltaT" : 10},
		{"year" : 1709 , "deltaT" : 10},
		{"year" : 1710 , "deltaT" : 10},
		{"year" : 1711 , "deltaT" : 10},
		{"year" : 1712 , "deltaT" : 10},
		{"year" : 1713 , "deltaT" : 10},
		{"year" : 1714 , "deltaT" : 10},
		{"year" : 1715 , "deltaT" : 10},
		{"year" : 1716 , "deltaT" : 10},
		{"year" : 1717 , "deltaT" : 11},
		{"year" : 1718 , "deltaT" : 11},
		{"year" : 1719 , "deltaT" : 11},
		{"year" : 1720 , "deltaT" : 11},
		{"year" : 1721 , "deltaT" : 11},
		{"year" : 1722 , "deltaT" : 11},
		{"year" : 1723 , "deltaT" : 11},
		{"year" : 1724 , "deltaT" : 11},
		{"year" : 1725 , "deltaT" : 11},
		{"year" : 1726 , "deltaT" : 11},
		{"year" : 1727 , "deltaT" : 11},
		{"year" : 1728 , "deltaT" : 11},
		{"year" : 1729 , "deltaT" : 11},
		{"year" : 1730 , "deltaT" : 11},
		{"year" : 1731 , "deltaT" : 11},
		{"year" : 1732 , "deltaT" : 11},
		{"year" : 1733 , "deltaT" : 11},
		{"year" : 1734 , "deltaT" : 12},
		{"year" : 1735 , "deltaT" : 12},
		{"year" : 1736 , "deltaT" : 12},
		{"year" : 1737 , "deltaT" : 12},
		{"year" : 1738 , "deltaT" : 12},
		{"year" : 1739 , "deltaT" : 12},
		{"year" : 1740 , "deltaT" : 12},
		{"year" : 1741 , "deltaT" : 12},
		{"year" : 1742 , "deltaT" : 12},
		{"year" : 1743 , "deltaT" : 12},
		{"year" : 1744 , "deltaT" : 13},
		{"year" : 1745 , "deltaT" : 13},
		{"year" : 1746 , "deltaT" : 13},
		{"year" : 1747 , "deltaT" : 13},
		{"year" : 1748 , "deltaT" : 13},
		{"year" : 1749 , "deltaT" : 13},
		{"year" : 1750 , "deltaT" : 13},
		{"year" : 1751 , "deltaT" : 14},
		{"year" : 1752 , "deltaT" : 14},
		{"year" : 1753 , "deltaT" : 14},
		{"year" : 1754 , "deltaT" : 14},
		{"year" : 1755 , "deltaT" : 14},
		{"year" : 1756 , "deltaT" : 14},
		{"year" : 1757 , "deltaT" : 14},
		{"year" : 1758 , "deltaT" : 15},
		{"year" : 1759 , "deltaT" : 15},
		{"year" : 1760 , "deltaT" : 15},
		{"year" : 1761 , "deltaT" : 15},
		{"year" : 1762 , "deltaT" : 15},
		{"year" : 1763 , "deltaT" : 15},
		{"year" : 1764 , "deltaT" : 15},
		{"year" : 1765 , "deltaT" : 16},
		{"year" : 1766 , "deltaT" : 16},
		{"year" : 1767 , "deltaT" : 16},
		{"year" : 1768 , "deltaT" : 16},
		{"year" : 1769 , "deltaT" : 16},
		{"year" : 1770 , "deltaT" : 16},
		{"year" : 1771 , "deltaT" : 16},
		{"year" : 1772 , "deltaT" : 16},
		{"year" : 1773 , "deltaT" : 16},
		{"year" : 1774 , "deltaT" : 16},
		{"year" : 1775 , "deltaT" : 17},
		{"year" : 1776 , "deltaT" : 17},
		{"year" : 1777 , "deltaT" : 17},
		{"year" : 1778 , "deltaT" : 17},
		{"year" : 1779 , "deltaT" : 17},
		{"year" : 1780 , "deltaT" : 17},
		{"year" : 1781 , "deltaT" : 17},
		{"year" : 1782 , "deltaT" : 17},
		{"year" : 1783 , "deltaT" : 17},
		{"year" : 1784 , "deltaT" : 17},
		{"year" : 1785 , "deltaT" : 17},
		{"year" : 1786 , "deltaT" : 17},
		{"year" : 1787 , "deltaT" : 17},
		{"year" : 1788 , "deltaT" : 17},
		{"year" : 1789 , "deltaT" : 17},
		{"year" : 1790 , "deltaT" : 17},
		{"year" : 1791 , "deltaT" : 17},
		{"year" : 1792 , "deltaT" : 16},
		{"year" : 1793 , "deltaT" : 16},
		{"year" : 1794 , "deltaT" : 16},
		{"year" : 1795 , "deltaT" : 16},
		{"year" : 1796 , "deltaT" : 15},
		{"year" : 1797 , "deltaT" : 15},
		{"year" : 1798 , "deltaT" : 14},
		{"year" : 1799 , "deltaT" : 14},
		{"year" : 1800 , "deltaT" : 13.7},
		{"year" : 1801 , "deltaT" : 13.4},
		{"year" : 1802 , "deltaT" : 13.1},
		{"year" : 1803 , "deltaT" : 12.9},
		{"year" : 1804 , "deltaT" : 12.7},
		{"year" : 1805 , "deltaT" : 12.6},
		{"year" : 1806 , "deltaT" : 12.5},
		{"year" : 1807 , "deltaT" : 12.5},
		{"year" : 1808 , "deltaT" : 12.5},
		{"year" : 1809 , "deltaT" : 12.5},
		{"year" : 1810 , "deltaT" : 12.5},
		{"year" : 1811 , "deltaT" : 12.5},
		{"year" : 1812 , "deltaT" : 12.5},
		{"year" : 1813 , "deltaT" : 12.5},
		{"year" : 1814 , "deltaT" : 12.5},
		{"year" : 1815 , "deltaT" : 12.5},
		{"year" : 1816 , "deltaT" : 12.5},
		{"year" : 1817 , "deltaT" : 12.4},
		{"year" : 1818 , "deltaT" : 12.3},
		{"year" : 1819 , "deltaT" : 12.2},
		{"year" : 1820 , "deltaT" : 12.0},
		{"year" : 1821 , "deltaT" : 11.7},
		{"year" : 1822 , "deltaT" : 11.4},
		{"year" : 1823 , "deltaT" : 11.1},
		{"year" : 1824 , "deltaT" : 10.6},
		{"year" : 1825 , "deltaT" : 10.2},
		{"year" : 1826 , "deltaT" : 9.6},
		{"year" : 1827 , "deltaT" : 9.1},
		{"year" : 1828 , "deltaT" : 8.6},
		{"year" : 1829 , "deltaT" : 8.0},
		{"year" : 1830 , "deltaT" : 7.5},
		{"year" : 1831 , "deltaT" : 7.0},
		{"year" : 1832 , "deltaT" : 6.6},
		{"year" : 1833 , "deltaT" : 6.3},
		{"year" : 1834 , "deltaT" : 6.0},
		{"year" : 1835 , "deltaT" : 5.8},
		{"year" : 1836 , "deltaT" : 5.7},
		{"year" : 1837 , "deltaT" : 5.6},
		{"year" : 1838 , "deltaT" : 5.6},
		{"year" : 1839 , "deltaT" : 5.6},
		{"year" : 1840 , "deltaT" : 5.7},
		{"year" : 1841 , "deltaT" : 5.8},
		{"year" : 1842 , "deltaT" : 5.9},
		{"year" : 1843 , "deltaT" : 6.1},
		{"year" : 1844 , "deltaT" : 6.2},
		{"year" : 1845 , "deltaT" : 6.3},
		{"year" : 1846 , "deltaT" : 6.5},
		{"year" : 1847 , "deltaT" : 6.6},
		{"year" : 1848 , "deltaT" : 6.8},
		{"year" : 1849 , "deltaT" : 6.9},
		{"year" : 1850 , "deltaT" : 7.1},
		{"year" : 1851 , "deltaT" : 7.2},
		{"year" : 1852 , "deltaT" : 7.3},
		{"year" : 1853 , "deltaT" : 7.4},
		{"year" : 1854 , "deltaT" : 7.5},
		{"year" : 1855 , "deltaT" : 7.6},
		{"year" : 1856 , "deltaT" : 7.7},
		{"year" : 1857 , "deltaT" : 7.7},
		{"year" : 1858 , "deltaT" : 7.8},
		{"year" : 1859 , "deltaT" : 7.8},
		{"year" : 1860 , "deltaT" : 7.88},
		{"year" : 1861 , "deltaT" : 7.82},
		{"year" : 1862 , "deltaT" : 7.54},
		{"year" : 1863 , "deltaT" : 6.97},
		{"year" : 1864 , "deltaT" : 6.40},
		{"year" : 1865 , "deltaT" : 6.02},
		{"year" : 1866 , "deltaT" : 5.41},
		{"year" : 1867 , "deltaT" : 4.10},
		{"year" : 1868 , "deltaT" : 2.92},
		{"year" : 1869 , "deltaT" : 1.82},
		{"year" : 1870 , "deltaT" : 1.61},
		{"year" : 1871 , "deltaT" : 0.10},
		{"year" : 1872 , "deltaT" : -1.02},
		{"year" : 1873 , "deltaT" : -1.28},
		{"year" : 1874 , "deltaT" : -2.69},
		{"year" : 1875 , "deltaT" : -3.24},
		{"year" : 1876 , "deltaT" : -3.64},
		{"year" : 1877 , "deltaT" : -4.54},
		{"year" : 1878 , "deltaT" : -4.71},
		{"year" : 1879 , "deltaT" : -5.11},
		{"year" : 1880 , "deltaT" : -5.40},
		{"year" : 1881 , "deltaT" : -5.42},
		{"year" : 1882 , "deltaT" : -5.20},
		{"year" : 1883 , "deltaT" : -5.46},
		{"year" : 1884 , "deltaT" : -5.46},
		{"year" : 1885 , "deltaT" : -5.79},
		{"year" : 1886 , "deltaT" : -5.63},
		{"year" : 1887 , "deltaT" : -5.64},
		{"year" : 1888 , "deltaT" : -5.80},
		{"year" : 1889 , "deltaT" : -5.66},
		{"year" : 1890 , "deltaT" : -5.87},
		{"year" : 1891 , "deltaT" : -6.01},
		{"year" : 1892 , "deltaT" : -6.19},
		{"year" : 1893 , "deltaT" : -6.64},
		{"year" : 1894 , "deltaT" : -6.44},
		{"year" : 1895 , "deltaT" : -6.47},
		{"year" : 1896 , "deltaT" : -6.09},
		{"year" : 1897 , "deltaT" : -5.76},
		{"year" : 1898 , "deltaT" : -4.66},
		{"year" : 1899 , "deltaT" : -3.74},
		{"year" : 1900 , "deltaT" : -2.72},
		{"year" : 1901 , "deltaT" : -1.54},
		{"year" : 1902 , "deltaT" : -0.02},
		{"year" : 1903 , "deltaT" : 1.24},
		{"year" : 1904 , "deltaT" : 2.64},
		{"year" : 1905 , "deltaT" : 3.86},
		{"year" : 1906 , "deltaT" : 5.37},
		{"year" : 1907 , "deltaT" : 6.14},
		{"year" : 1908 , "deltaT" : 7.75},
		{"year" : 1909 , "deltaT" : 9.13},
		{"year" : 1910 , "deltaT" : 10.46},
		{"year" : 1911 , "deltaT" : 11.53},
		{"year" : 1912 , "deltaT" : 13.36},
		{"year" : 1913 , "deltaT" : 14.65},
		{"year" : 1914 , "deltaT" : 16.01},
		{"year" : 1915 , "deltaT" : 17.20},
		{"year" : 1916 , "deltaT" : 18.24},
		{"year" : 1917 , "deltaT" : 19.06},
		{"year" : 1918 , "deltaT" : 20.25},
		{"year" : 1919 , "deltaT" : 20.95},
		{"year" : 1920 , "deltaT" : 21.16},
		{"year" : 1921 , "deltaT" : 22.25},
		{"year" : 1922 , "deltaT" : 22.41},
		{"year" : 1923 , "deltaT" : 23.03},
		{"year" : 1924 , "deltaT" : 23.49},
		{"year" : 1925 , "deltaT" : 23.62},
		{"year" : 1926 , "deltaT" : 23.68},
		{"year" : 1927 , "deltaT" : 24.49},
		{"year" : 1928 , "deltaT" : 24.34},
		{"year" : 1929 , "deltaT" : 24.08},
		{"year" : 1930 , "deltaT" : 24.02},
		{"year" : 1931 , "deltaT" : 24.00},
		{"year" : 1932 , "deltaT" : 23.87},
		{"year" : 1933 , "deltaT" : 23.95},
		{"year" : 1934 , "deltaT" : 23.86},
		{"year" : 1935 , "deltaT" : 23.93},
		{"year" : 1936 , "deltaT" : 23.73},
		{"year" : 1937 , "deltaT" : 23.92},
		{"year" : 1938 , "deltaT" : 23.96},
		{"year" : 1939 , "deltaT" : 24.02},
		{"year" : 1940 , "deltaT" : 24.33},
		{"year" : 1941 , "deltaT" : 24.83},
		{"year" : 1942 , "deltaT" : 25.30},
		{"year" : 1943 , "deltaT" : 25.70},
		{"year" : 1944 , "deltaT" : 26.24},
		{"year" : 1945 , "deltaT" : 26.77},
		{"year" : 1946 , "deltaT" : 27.28},
		{"year" : 1947 , "deltaT" : 27.78},
		{"year" : 1948 , "deltaT" : 28.25},
		{"year" : 1949 , "deltaT" : 28.71},
		{"year" : 1950 , "deltaT" : 29.15},
		{"year" : 1951 , "deltaT" : 29.57},
		{"year" : 1952 , "deltaT" : 29.97},
		{"year" : 1953 , "deltaT" : 30.36},
		{"year" : 1954 , "deltaT" : 30.72},
		{"year" : 1955 , "deltaT" : 31.07},
		{"year" : 1956 , "deltaT" : 31.35},
		{"year" : 1957 , "deltaT" : 31.68},
		{"year" : 1958 , "deltaT" : 32.18},
		{"year" : 1959 , "deltaT" : 32.68},
		{"year" : 1960 , "deltaT" : 33.15},
		{"year" : 1961 , "deltaT" : 33.59},
		{"year" : 1962 , "deltaT" : 34.00},
		{"year" : 1963 , "deltaT" : 34.47},
		{"year" : 1964 , "deltaT" : 35.03},
		{"year" : 1965 , "deltaT" : 35.73},
		{"year" : 1966 , "deltaT" : 36.54},
		{"year" : 1967 , "deltaT" : 37.43},
		{"year" : 1968 , "deltaT" : 38.29},
		{"year" : 1969 , "deltaT" : 39.20},
		{"year" : 1970 , "deltaT" : 40.18},
		{"year" : 1971 , "deltaT" : 41.17},
		{"year" : 1972 , "deltaT" : 42.23},
		{"year" : 1973 , "deltaT" : 43.37},
		{"year" : 1974 , "deltaT" : 44.4841},
		{"year" : 1975 , "deltaT" : 45.4761},
		{"year" : 1976 , "deltaT" : 46.4567},
		{"year" : 1977 , "deltaT" : 47.5214},
		{"year" : 1978 , "deltaT" : 48.5344},
		{"year" : 1979 , "deltaT" : 49.5861},
		{"year" : 1980 , "deltaT" : 50.5387},
		{"year" : 1981 , "deltaT" : 51.3808},
		{"year" : 1982 , "deltaT" : 52.1668},
		{"year" : 1983 , "deltaT" : 52.9565},
		{"year" : 1984 , "deltaT" : 53.7882},
		{"year" : 1985 , "deltaT" : 54.3427},
		{"year" : 1986 , "deltaT" : 54.8712},
		{"year" : 1987 , "deltaT" : 55.3222},
		{"year" : 1988 , "deltaT" : 55.8197},
		{"year" : 1989 , "deltaT" : 56.3000},
		{"year" : 1990 , "deltaT" : 56.8553},
		{"year" : 1991 , "deltaT" : 57.5653},
		{"year" : 1992 , "deltaT" : 58.3092},
		{"year" : 1993 , "deltaT" : 59.1218},
		{"year" : 1994 , "deltaT" : 59.9845},
		{"year" : 1995 , "deltaT" : 60.7853},
		{"year" : 1996 , "deltaT" : 61.6287},
		{"year" : 1997 , "deltaT" : 62.2950},
		{"year" : 1998 , "deltaT" : 62.9659},
		{"year" : 1999 , "deltaT" : 63.4673},
		{"year" : 2000 , "deltaT" : 63.8285},
		{"year" : 2001 , "deltaT" : 64.0908},
		{"year" : 2002 , "deltaT" : 64.2998},
		{"year" : 2003 , "deltaT" : 64.4734},
		{"year" : 2004 , "deltaT" : 64.5736},
		{"year" : 2005 , "deltaT" : 64.6876},
		{"year" : 2006 , "deltaT" : 64.8452},
		{"year" : 2007 , "deltaT" : 65.1464},
		{"year" : 2008 , "deltaT" : 65.4574},
		{"year" : 2009 , "deltaT" : 65.7768},
		{"year" : 2010 , "deltaT" : 66.0699},
		{"year" : 2011 , "deltaT" : 66.3246},
		{"year" : 2012 , "deltaT" : 66.6030},
		{"year" : 2013 , "deltaT" : 66.9069},
		{"year" : 2014 , "deltaT" : 67.2810},
		{"year" : 2015 , "deltaT" : 67.6439},
		{"year" : 2016 , "deltaT" : 68.1024},
		{"year" : 2017 , "deltaT" : 68.5927},
		{"year" : 2018 , "deltaT" : 68.9677},
		{"year" : 2019 , "deltaT" : 69.2202},
		
		{"year" : 2020 , "deltaT" : 69.87},
		{"year" : 2021 , "deltaT" : 70.39},
		{"year" : 2022 , "deltaT" : 70.91},
		{"year" : 2023 , "deltaT" : 71.40},
		{"year" : 2024 , "deltaT" : 71.88},
		{"year" : 2025 , "deltaT" : 72.36},
		{"year" : 2026 , "deltaT" : 72.83},
		{"year" : 2027 , "deltaT" : 73.32},
	];
	let	diff = 0,
		bestIndex = 0,
		bestDiff = Infinity,
		i, cur;

	for (i = 0; i < deltaTable.length; i++) {
		cur = deltaTable[i];
		diff = Math.abs(year - cur.year);
		if (diff < bestDiff) {
			bestDiff = diff;
			bestIndex = i;
		}
	}
	
	x1 = bestIndex - 1;
	x2 = bestIndex;
	x3 = bestIndex + 1;
	
	y1 = deltaTable[x1].deltaT;
	y2 = deltaTable[x2].deltaT;
	y3 = deltaTable[x3].deltaT;
	
	a = y2 - y1;
	b = y3 - y2;
	c = b - a;
	n = bestDiff;
	
	let deltaT = y2 + n/2 * (a + b + n*c);
	
	let JDE = JD + deltaT/86400.0;
	return JDE;
}

function sunEquatorCoords(JD) {
	
}

//UTC date time
let year = 2012;
let month = 10;
let day = 12;
let hour = 0;
let min = 0;
let sec = 0;

let JD = julianDay(year, month, day, hour, min, sec);
let coords = [313.366800,-23.566600,0.0000000]; //Sao Paulo

//sunPosition(JD, coords);
sunPosition(JD, coords);

//let sun = sunHorizontPosition(JD, coordinates);

//let sun = sunEquatorCoords(JD);

//let date = JD2Date(JD, "year");
//let time = dayDec2Hour(date.day);


//document.write('Date = '+ date.year +'-'+ date.month +'-'+ date.day);
document.write('<BR>');
//document.write('time = '+ time.hour +':'+ time.min +':'+ time.sec);
document.write('<BR>');	

document.write('Julian Date = '+ JD);
document.write('<BR>');	
document.write('Date = '+ year +'-'+ month +'-'+ day +' '+ hour.toString().padStart(2,'0') +':'+ min.toString().padStart(2,'0') +':'+ sec.toString().padStart(2,'0'));
document.write('<BR>');	
document.write('<BR>');
//document.write('sun RA = '+ sun.RA);
document.write('<BR>');
//document.write('sun Declinacao = '+ sun.Declination);
document.write('<BR>');
document.write('<BR>');
//document.write('sun Azimuth = '+ sun.Azimuth);
document.write('<BR>');
//document.write('sun Altitude = '+ sun.Altitude);
document.write('<BR>');
//document.write('sun Hour Angle = '+ sun.HourAngle);
document.write('<BR>');
//document.write('sun LST = '+ sun.LST);
document.write('<BR>');