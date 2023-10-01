$fn = 96;

/* Library for involute gears, screwn and linear_rack

Contains the Module
- linear_rack(modul, length, height, width, contact_angle=20, helix_angle=0)
- spur_gear(modul, teeth, width, bore_hole, contact_angle=20, helix_angle=0, optimized=true)
- arrow_gear(modul, teeth, width, bore_hole, contact_angle=20, helix_angle=0, optimized=true)
- rack_and_wheel (modul, length_stange, teeth_rad, height_stange, bore_hole_rad, width, contact_angle=20, helix_angle=0, assembled=true, optimized=true)
- angled_inner_gear(modul, teeth, width, randwidth, contact_angle=20, helix_angle=0)
- arrow_inner_gear(modul, teeth, width, randwidth, contact_angle=20, helix_angle=0)
- planetary_gear(modul, teeth_sonne, teeth_planet, anzahl_planeten, width, randwidth, bore_hole, contact_angle=20, helix_angle=0, assembled=true, optimized=true)
- bevel_gear(modul, teeth,  cone_angle, gear_width, bore_hole, contact_angle=20, helix_angle=0)
- arrow_bevel_gear(modul, teeth, cone_angle, gear_width, bore_hole, contact_angle=20, helix_angle=0)
- cone_wheel_pair(modul, teeth_rad, teeth_pinion, axis_angle=90, gear_width, bore_hole, contact_angle = 20, helix_angle=0, assembled=true)
- locking_wheel_pair(modul, teeth_rad, number_of_teeth, axis_angle=90, gear_width, bore_hole, contact_angle = 20, helix_angle=0, assembled=true)
- screw(modul, gear_number, length, bore_hole, contact_angle=20, pitch_angle=10, assembled=true)
- screw_gear_set(modul, teeth, gear_number, width, length, bore_hole, bore_hole_rad, contact_angle=20, pitch_angle=0, optimized=true, assembled=true)

Examples for each module are commented out at the end of this file

Autor:		Dr Jörg Janssen
Stand:		29. Oktober 2018
Version:	2.3
Lizenz:		Creative Commons - Attribution, Non Commercial, Share Alike

Permitted modules according to DIN 780:
0.05 0.06 0.08 0.10 0.12 0.16
0.20 0.25 0.3  0.4  0.5  0.6
0.7  0.8  0.9  1    1.25 1.5
2    2.5  3    4    5    6
8    10   12   16   20   25
32   40   50   60

*/


// Common variables
pi = 3.14159;
rad = 57.29578;
spiel = 0.05;	// Play between teeth

/*	Converts radian to degrees */
function grad(contact_angle) = contact_angle*rad;

/*	Converts degrees to radian */
function radian(contact_angle) = contact_angle/rad;

/*	Converts 2D polar coordinates into Cartesian ones
    Format: radius, phi; phi = angle to x-axis on xy-plane */
function pole2point(polvect) = [
	polvect[0]*cos(polvect[1]),  
	polvect[0]*sin(polvect[1])
];

/*	Circular involute function:
    Outputs the polar coordinates of a circular involute.
    r = Radius of the base circle
    rho = roll angle in degrees */
function ev(r,rho) = [
	r/cos(rho),
	grad(tan(rho)-radian(rho))
];

/*  Spherical involute function
    Outputs the azimuth angle of a spherical involute
    theta0 = polar angle of the cone at whose intersection edge to the large sphere the involute rolls off
    theta = polar angle for which the azimuth angle of the involute is to be calculated */
function sphere_calc(theta0,theta) = 1/sin(theta0)*acos(cos(theta)/cos(theta0))-acos(tan(theta0)/tan(theta));

/*  Converts spherical coordinates to Cartesian
    Format: radius, theta, phi; theta = angle to z-axis, phi = angle to x-axis on xy plane */
function sphere_xlate(vect) = [
	vect[0]*sin(vect[1])*cos(vect[2]),  
	vect[0]*sin(vect[1])*sin(vect[2]),
	vect[0]*cos(vect[1])
];

/*	checks if a number is even
	= 1, if yes
	= 0, if the number is not even */
function is_straight(zahl) =
	(zahl == floor(zahl/2)*2) ? 1 : 0;

/*	greatest common divisor
	according to Euclidean algorithm.
	Sorting: a must be greater than b */
function ggt(a,b) = 
	a%b == 0 ? b : ggt(b,a%b);

/*	Polar function with polar angle and two variables */
function spirale(a, r0, phi) =
	a*phi + r0; 

/*	Copies and rotates a body */
module copy_rot(vect, zahl, abstand, winkel){
	for(i = [0:zahl-1]){
		translate(v=vect*abstand*i)
			rotate(a=i*winkel, v = [0,0,1])
				children(0);
	}
}

/*  linear_rack
    modul = height of the tooth tip above the rolling line
    length = length of the linear_rack
    height = height of the linear_rack up to the pitch line
    width = width of a tooth
    contact_angle = contact_angle, default value = 20° according to DIN 867. Should not be greater than 45°.
    helix_angle = helix angle to the linear_rackn transverse axis; 0° = spur toothing */
module linear_rack(modul, length, height, width, contact_angle = 20, helix_angle = 0) {

	// Dimension calculations
	modul=modul*(1-spiel);
	c = modul / 6;												// Header
	mx = modul/cos(helix_angle);							// Modulus distorted by helix angle in x-direction
	a = 2*mx*tan(contact_angle)+c*tan(contact_angle);		// Edge width
	b = pi*mx/2-2*mx*tan(contact_angle);						// Header width
	x = width*tan(helix_angle);							// Displacement of the top side in x-direction due to helix angle
	nz = ceil((length+abs(2*x))/(pi*mx));						// Number of teeth
	
	translate([-pi*mx*(nz-1)/2-a-b/2,-modul,0]){
		intersection(){											// Creates a prism that is fitted into a cuboid geometry
			copy_rot([1,0,0], nz, pi*mx, 0){
				polyhedron(
					points=[[0,-c,0], [a,2*modul,0], [a+b,2*modul,0], [2*a+b,-c,0], [pi*mx,-c,0], [pi*mx,modul-height,0], [0,modul-height,0],	// Bottom
						[0+x,-c,width], [a+x,2*modul,width], [a+b+x,2*modul,width], [2*a+b+x,-c,width], [pi*mx+x,-c,width], [pi*mx+x,modul-height,width], [0+x,modul-height,width]],	// Top
					faces=[[6,5,4,3,2,1,0],						// Bottom
						[1,8,7,0],
						[9,8,1,2],
						[10,9,2,3],
						[11,10,3,4],
						[12,11,4,5],
						[13,12,5,6],
						[7,13,6,0],
						[7,8,9,10,11,12,13],					// Top
					]
				);
			};
			translate([abs(x),-height+modul-0.5,-0.5]){
				cube([length,height+modul+1,width+1]);			// Cuboid comprising the volume of the linear_rack
			}	
		};
	};	
}

/*  spur_gear
    modul = height of the tooth tip above the pitch circle
    teeth = number of gear teeth
    width = gear_width
    bore_hole = diameter of the center bore_hole
    contact_angle = contact_angle, default value = 20° according to DIN 867. should not be greater than 45°.
    helix_angle = helix angle to the axis of rotation; 0° = spur toothing
	optimized = Holes to save material/weight or create surface enlargement, if geometry allows. */
module spur_gear(modul, teeth, width, bore_hole, contact_angle = 20, helix_angle = 0, optimized = true) {

	// Dimension calculations	
	d = modul * teeth;											// Pitch circle diameter
	r = d / 2;														// Pitch radius
	alpha_stirn = atan(tan(contact_angle)/cos(helix_angle));// Bevel angle in face cut
	db = d * cos(alpha_stirn);										// Base circle diameter
	rb = db / 2;													// Base circle radius
	da = (modul <1)? d + modul * 2.2 : d + modul * 2;				// tip diameter according to DIN 58400 or DIN 867
	ra = da / 2;													// Tip circle radius
	c =  (teeth <3)? 0 : modul/6;								// Header play
	df = d - 2 * (modul + c);										// Base circle diameter
	rf = df / 2;													// Foot circle radius
	rho_ra = acos(rb/ra);											// maximum rolling angle;
																	// Involute starts at base circle and ends at tip circle
	rho_r = acos(rb/r);												// Roll angle on the pitch circle;
																	// Roll angle on the pitch circle;
	phi_r = grad(tan(rho_r)-radian(rho_r));							// Angle to the point of the involute on pitch circle
	gamma = rad*width/(r*tan(90-helix_angle));				// Torsion angle for extrusion
	a_step = rho_ra/16;											// Involute is divided into 16 pieces
	tau = 360/teeth;												// Angle of pitch
	
	r_loch = (2*rf - bore_hole)/8;									// Radius of the holes für Material-/Gewichtsersparnis
	rm = bore_hole/2+2*r_loch;										// Distance of the axes of the holes from the main axis
	z_loch = floor(2*pi*rm/(3*r_loch));								// Number of holes for material/weight saving
	
	optimized = (optimized && r >= width*1.5 && d > 2*bore_hole);	// is optimization useful?

	// Drawing
	union(){
		rotate([0,0,-phi_r-90*(1-spiel)/teeth]){						// Center tooth on x-axis;
																		// Makes alignment with other wheels easier

			linear_extrude(height = width, twist = gamma){
				difference(){
					union(){
						gear_width = (180*(1-spiel))/teeth+2*phi_r;
						circle(rf);										// Fußkreis	
						for (rot = [0:tau:360]){
							rotate (rot){								// "teeth-mal" copy_rotn und drehen
								polygon(concat(							// Zahn
									[[0,0]],							// Tooth segment begins and ends at the origin
									[for (rho = [0:a_step:rho_ra])		// von null Grad (Base Circle)
																		// up to maximum involute angle (tip circle)
										pole2point(ev(rb,rho))],		// First involute edge

									[pole2point(ev(rb,rho_ra))],		// Point of involute on tip circle

									[for (rho = [rho_ra:-a_step:0])	// of maximum involute angle (tip circle)
																		// to zero degrees (Base Circle)
										pole2point([ev(rb,rho)[0], gear_width-ev(rb,rho)[1]])]
																		// Second involute edge
																		// (180*(1-play)) instead of 180 degrees,
																		// to allow play on the flanks
									)
								);
							}
						}
					}			
					circle(r = rm+r_loch*1.49);							// "bore_hole"
				}
			}
		}
		// with material savings
		if (optimized) {
			linear_extrude(height = width){
				difference(){
						circle(r = (bore_hole+r_loch)/2);
						circle(r = bore_hole/2);							// bore_hole
					}
				}
			linear_extrude(height = (width-r_loch/2 < width*2/3) ? width*2/3 : width-r_loch/2){
				difference(){
					circle(r=rm+r_loch*1.51);
					union(){
						circle(r=(bore_hole+r_loch)/2);
						for (i = [0:1:z_loch]){
							translate(sphere_xlate([rm,90,i*360/z_loch]))
								circle(r = r_loch);
						}
					}
				}
			}
		}
		// without saving material
		else {
			linear_extrude(height = width){
				difference(){
					circle(r = rm+r_loch*1.51);
					circle(r = bore_hole/2);
				}
			}
		}
	}
}

/*  arrow_gear; uses the module "spur_gear"
    modul = Height of the tooth head above the pitch circle
    teeth = Number of wheel teeth
    width = gear_width
    bore_hole = Center hole diameter
    contact_angle = contact_angle, Standard value = 20° according to DIN 867. Should not be greater than 45° sein.
    helix_angle = Helix angle to the axis of rotation, standard value = 0° (straight teeth)
	optimized = Holes to save material/weight */
module arrow_gear(modul, teeth, width, bore_hole, contact_angle = 20, helix_angle=0, optimized=true){

	width = width/2;
	d = modul * teeth;											// Pitch circle diameter
	r = d / 2;														// Pitch radius
	c =  (teeth <3)? 0 : modul/6;								// Header play

	df = d - 2 * (modul + c);										// Base circle diameter
	rf = df / 2;													// Foot circle radius

	r_loch = (2*rf - bore_hole)/8;									// Radius of the holes für Material-/Gewichtsersparnis
	rm = bore_hole/2+2*r_loch;										// Distance of the axes of the holes from the main axis
	z_loch = floor(2*pi*rm/(3*r_loch));								// Number of holes for material/weight savings
	
	optimized = (optimized && r >= width*3 && d > 2*bore_hole);		// is optimization useful?

	translate([0,0,width]){
		union(){
			spur_gear(modul, teeth, width, 2*(rm+r_loch*1.49), contact_angle, helix_angle, false);		// untere Hälfte
			mirror([0,0,1]){
				spur_gear(modul, teeth, width, 2*(rm+r_loch*1.49), contact_angle, helix_angle, false);	// obere Hälfte
			}
		}
	}
	// with material savings
	if (optimized) {
		linear_extrude(height = width*2){
			difference(){
					circle(r = (bore_hole+r_loch)/2);
					circle(r = bore_hole/2);							// bore_hole
				}
			}
		linear_extrude(height = (2*width-r_loch/2 < 1.33*width) ? 1.33*width : 2*width-r_loch/2){ //width*4/3
			difference(){
				circle(r=rm+r_loch*1.51);
				union(){
					circle(r=(bore_hole+r_loch)/2);
					for (i = [0:1:z_loch]){
						translate(sphere_xlate([rm,90,i*360/z_loch]))
							circle(r = r_loch);
					}
				}
			}
		}
	}
	// without saving material
	else {
		linear_extrude(height = width*2){
			difference(){
				circle(r = rm+r_loch*1.51);
				circle(r = bore_hole/2);
			}
		}
	}
}

/*	linear_rack und -Rad
    modul = Height of the tooth head above the pitch circle
    length_stange = length der linear_rack
    teeth_rad = Number of wheel teeth
	height_stange = Height of linear_rack bis zur Wälzgeraden
    bore_hole_rad = Center hole diameter des spur_gears
	width = width eines Zahns
    contact_angle = contact_angle, Standard value = 20° according to DIN 867. Should not be greater than 45° sein.
    helix_angle = Helix angle to the axis of rotation, standard value = 0° (straight teeth) */
module rack_and_wheel (modul, length_stange, teeth_rad, height_stange, bore_hole_rad, width, contact_angle=20, helix_angle=0, assembled=true, optimized=true) {

	abstand = assembled? modul*teeth_rad/2 : modul*teeth_rad;

	linear_rack(modul, length_stange, height_stange, width, contact_angle, -helix_angle);
	translate([0,abstand,0])
		rotate(a=360/teeth_rad)
			spur_gear (modul, teeth_rad, width, bore_hole_rad, contact_angle, helix_angle, optimized);
}

/*	angled_inner_gear
    modul = Height of the tooth head above the pitch circle
    teeth = Number of wheel teeth
    width = gear_width
	randwidth = width des Randes ab Fußkreis
    bore_hole = Center hole diameter
    contact_angle = contact_angle, Standard value = 20° according to DIN 867. Should not be greater than 45° sein.
    helix_angle = Helix angle to the axis of rotation, standard value = 0° (straight teeth) */
module angled_inner_gear(modul, teeth, width, randwidth, contact_angle = 20, helix_angle = 0) {

	// Dimension calculations	
	ha = (teeth >= 20) ? 0.02 * atan((teeth/15)/pi) : 0.6;	// Shortening factor of tooth head height
	d = modul * teeth;											// Pitch circle diameter
	r = d / 2;														// Pitch radius
	alpha_stirn = atan(tan(contact_angle)/cos(helix_angle));// Bevel angle in face cut
	db = d * cos(alpha_stirn);										// Base circle diameter
	rb = db / 2;													// Base circle radius
	c = modul / 6;													// Header play
	da = (modul <1)? d + (modul+c) * 2.2 : d + (modul+c) * 2;		// Kopfkreisdurchmesser
	ra = da / 2;													// Tip circle radius
	df = d - 2 * modul * ha;										// Base circle diameter
	rf = df / 2;													// Foot circle radius
	rho_ra = acos(rb/ra);											// maximaler Evolventenwinkel;
																	// Evolvente beginnt auf Base Circle und endet an Kopfkreis
	rho_r = acos(rb/r);												// Evolventenwinkel am Teilkreis;
																	// Evolvente beginnt auf Base Circle und endet an Kopfkreis
	phi_r = grad(tan(rho_r)-radian(rho_r));							// Winkel zum Punkt der Evolvente auf Teilkreis
	gamma = rad*width/(r*tan(90-helix_angle));				// Torsionswinkel für Extrusion
	a_step = rho_ra/16;											// Evolvente wird in 16 Stücke geteilt
	tau = 360/teeth;												// Angle of pitch

	// Drawing
	rotate([0,0,-phi_r-90*(1+spiel)/teeth])						// Center tooth on x-axis;
																	// Makes alignment with other wheels easier
	linear_extrude(height = width, twist = gamma){
		difference(){
			circle(r = ra + randwidth);							// Außenkreis
			union(){
				gear_width = (180*(1+spiel))/teeth+2*phi_r;
				circle(rf);											// Fußkreis	
				for (rot = [0:tau:360]){
					rotate (rot) {									// "teeth-mal" copy_rotn und drehen
						polygon( concat(
							[[0,0]],
							[for (rho = [0:a_step:rho_ra])			// von null Grad (Base Circle)
																	// bis maximaler Evolventenwinkel (Kopfkreis)
								pole2point(ev(rb,rho))],
							[pole2point(ev(rb,rho_ra))],
							[for (rho = [rho_ra:-a_step:0])		// von maximaler Evolventenwinkel (Kopfkreis)
																	// to zero degrees (Base Circle)
								pole2point([ev(rb,rho)[0], gear_width-ev(rb,rho)[1]])]
																	// (180*(1+spiel)) statt 180,
																	// to allow play on the flanks
							)
						);
					}
				}
			}
		}
	}

	echo("Außendurchmesser angled_inner_gear = ", 2*(ra + randwidth));
	
}

/*  Pfeil-angled_inner_gear; uses the module "angled_inner_gear"
    modul = Höhe des Zahnkopfes über dem Teilkegel
    teeth = Number of wheel teeth
    width = gear_width
    bore_hole = Center hole diameter
    contact_angle = contact_angle, Standard value = 20° according to DIN 867. Should not be greater than 45° sein.
    helix_angle = Helix angle to the axis of rotation, standard value = 0° (straight teeth) */
module arrow_inner_gear(modul, teeth, width, randwidth, contact_angle = 20, helix_angle = 0) {

	width = width / 2;
	translate([0,0,width])
		union(){
		angled_inner_gear(modul, teeth, width, randwidth, contact_angle, helix_angle);		// untere Hälfte
		mirror([0,0,1])
			angled_inner_gear(modul, teeth, width, randwidth, contact_angle, helix_angle);	// obere Hälfte
	}
}

/*	planetary_gear; verwendet die Module "arrow_gear" und "arrow_inner_gear"
    modul = Höhe des Zahnkopfes über dem Teilkegel
    teeth_sonne = Anzahl der Zähne des Sonnenrads
    teeth_planet = Anzahl der Zähne eines Planetenrads
    anzahl_planeten = Anzahl der Planetenräder. Wenn null, rechnet die Funktion die Mindestanzahl aus.
    width = gear_width
	randwidth = width des Randes ab Fußkreis
    bore_hole = Center hole diameter
    contact_angle = contact_angle, Standard value = 20° according to DIN 867. Should not be greater than 45° sein.
    helix_angle = Helix angle to the axis of rotation, standard value = 0° (straight teeth)
	assembled = 
	optimized = Holes to save material/weight bzw. Oberflächenvergößerung erzeugen, wenn Geometrie erlaubt
	assembled = Komponenten zusammengebaut für Konstruktion oder auseinander zum 3D-Druck	*/
module planetary_gear(modul, teeth_sonne, teeth_planet, anzahl_planeten, width, randwidth, bore_hole, contact_angle=20, helix_angle=0, assembled=true, optimized=true){

	// Dimension calculations
	d_sonne = modul*teeth_sonne;										// Pitch circle diameter Sonne
	d_planet = modul*teeth_planet;									// Pitch circle diameter Planeten
	achsabstand = modul*(teeth_sonne +  teeth_planet) / 2;		// Abstand von Sonnenrad-/angled_inner_gearachse und Planetenachse
	teeth_angled_inner_gear = teeth_sonne + 2*teeth_planet;				// Anzahl der Zähne des angled_inner_geares
    d_angled_inner_gear = modul*teeth_angled_inner_gear;									// Pitch circle diameter angled_inner_gear

	drehen = is_straight(teeth_planet);								// Muss das Sonnenrad gedreht werden?
		
	n_max = floor(180/asin(modul*(teeth_planet)/(modul*(teeth_sonne +  teeth_planet))));
																		// Anzahl Planetenräder: höchstens so viele, wie ohne
																		// Überlappung möglich

	// Drawing
	rotate([0,0,180/teeth_sonne*drehen]){
		arrow_gear (modul, teeth_sonne, width, bore_hole, contact_angle, -helix_angle, optimized);		// Sonnenrad
	}

	if (assembled){
        if(anzahl_planeten==0){
            list = [ for (n=[2 : 1 : n_max]) if ((((teeth_angled_inner_gear+teeth_sonne)/n)==floor((teeth_angled_inner_gear+teeth_sonne)/n))) n];
            anzahl_planeten = list[0];										// Ermittele Anzahl Planetenräder
             achsabstand = modul*(teeth_sonne + teeth_planet)/2;		// Abstand von Sonnenrad-/angled_inner_gearachse
            for(n=[0:1:anzahl_planeten-1]){
                translate(sphere_xlate([achsabstand,90,360/anzahl_planeten*n]))
					rotate([0,0,n*360*d_sonne/d_planet])
						arrow_gear (modul, teeth_planet, width, bore_hole, contact_angle, helix_angle);	// Planetenräder
            }
       }
       else{
            achsabstand = modul*(teeth_sonne + teeth_planet)/2;		// Abstand von Sonnenrad-/angled_inner_gearachse
            for(n=[0:1:anzahl_planeten-1]){
                translate(sphere_xlate([achsabstand,90,360/anzahl_planeten*n]))
                rotate([0,0,n*360*d_sonne/(d_planet)])
                    arrow_gear (modul, teeth_planet, width, bore_hole, contact_angle, helix_angle);	// Planetenräder
            }
		}
	}
	else{
		planetenabstand = teeth_angled_inner_gear*modul/2+randwidth+d_planet;		// Abstand Planeten untereinander
		for(i=[-(anzahl_planeten-1):2:(anzahl_planeten-1)]){
			translate([planetenabstand, d_planet*i,0])
				arrow_gear (modul, teeth_planet, width, bore_hole, contact_angle, helix_angle);	// Planetenräder
		}
	}

	arrow_inner_gear (modul, teeth_angled_inner_gear, width, randwidth, contact_angle, helix_angle); // angled_inner_gear

}

/*  bevel_gear
    modul = Höhe des Zahnkopfes über dem Teilkegel; Angabe für die Aussenseite des Kegels
    teeth = Number of wheel teeth
    cone_angle = (Halb)winkel des Kegels, auf dem das jeweils andere angled_inner_gear abrollt
    gear_width = width der Zähne von der Außenseite in Richtung Kegelspitze
    bore_hole = Center hole diameter
    contact_angle = contact_angle, Standard value = 20° according to DIN 867. Should not be greater than 45° sein.
	helix_angle = Schrägungswinkel, Standardwert = 0° */
module bevel_gear(modul, teeth, cone_angle, gear_width, bore_hole, contact_angle = 20, helix_angle=0) {

	// Dimension calculations
	d_outside = modul * teeth;									// Pitch cone diameter on the cone base surface,
																	// corresponds to the chord in the spherical section
	r_outside = d_outside / 2;										// Part of sphere radius auf der Kegelgrundfläche 
	rg_outside = r_outside/sin(cone_angle);						// Large sphere radius für Zahn-Außenseite, entspricht der Länge der Kegelflanke;
	rg_innen = rg_outside - gear_width;								// Large sphere radius für Zahn-Innenseite	
	r_innen = r_outside*rg_innen/rg_outside;
	alpha_stirn = atan(tan(contact_angle)/cos(helix_angle));// Bevel angle in face cut
	delta_b = asin(cos(alpha_stirn)*sin(cone_angle));			// Grundkegelwinkel		
	da_outside = (modul <1)? d_outside + (modul * 2.2) * cos(cone_angle): d_outside + modul * 2 * cos(cone_angle);
	ra_outside = da_outside / 2;
	delta_a = asin(ra_outside/rg_outside);
	c = modul / 6;													// Header play
	df_outside = d_outside - (modul +c) * 2 * cos(cone_angle);
	rf_outside = df_outside / 2;
	delta_f = asin(rf_outside/rg_outside);
	rkf = rg_outside*sin(delta_f);									// Radius des Kegelfußes
	height_f = rg_outside*cos(delta_f);								// Höhe des Kegels vom Fußkegel
	
	echo("Pitch cone diameter on the cone base surface = ", d_outside);
	
	// Größen für Komplementär-Kegelstumpf
	height_k = (rg_outside-gear_width)/cos(cone_angle);			// Höhe des Komplementärkegels für richtige Zahnlänge
	rk = (rg_outside-gear_width)/sin(cone_angle);				// Fußradius des Komplementärkegels
	rfk = rk*height_k*tan(delta_f)/(rk+height_k*tan(delta_f));		// Kopfradius des Zylinders für 
																	// Komplementär-Kegelstumpf
	height_fk = rk*height_k/(height_k*tan(delta_f)+rk);				// height des Komplementär-Kegelstumpfs

	echo("Höhe bevel_gear = ", height_f-height_fk);
	
	phi_r = sphere_calc(delta_b, cone_angle);						// Winkel zum Punkt der Evolvente auf Teilkegel
		
	// Torsionswinkel gamma aus Schrägungswinkel
	gamma_g = 2*atan(gear_width*tan(helix_angle)/(2*rg_outside-gear_width));
	gamma = 2*asin(rg_outside/r_outside*sin(gamma_g/2));
	
	a_step = (delta_a - delta_b)/16;
	tau = 360/teeth;												// Angle of pitch
	start = (delta_b > delta_f) ? delta_b : delta_f;
	spiegelpunkt = (180*(1-spiel))/teeth+2*phi_r;

	// Drawing
	rotate([0,0,phi_r+90*(1-spiel)/teeth]){						// Center tooth on x-axis;
																	// Makes alignment with other wheels easier
		translate([0,0,height_f]) rotate(a=[0,180,0]){
			union(){
				translate([0,0,height_f]) rotate(a=[0,180,0]){								// Kegelstumpf							
					difference(){
						linear_extrude(height=height_f-height_fk, scale=rfk/rkf) circle(rkf*1.001); // 1 promille Überlappung mit Zahnfuß
						translate([0,0,-1]){
							cylinder(h = height_f-height_fk+2, r = bore_hole/2);				// bore_hole
						}
					}	
				}
				for (rot = [0:tau:360]){
					rotate (rot) {															// "teeth-mal" copy_rotn und drehen
						union(){
							if (delta_b > delta_f){
								// Zahnfuß
								flankenpunkt_unten = 1*spiegelpunkt;
								flankenpunkt_oben = sphere_calc(delta_f, start);
								polyhedron(
									points = [
										sphere_xlate([rg_outside, start*1.001, flankenpunkt_unten]),	// 1 promille Überlappung mit Zahn
										sphere_xlate([rg_innen, start*1.001, flankenpunkt_unten+gamma]),
										sphere_xlate([rg_innen, start*1.001, spiegelpunkt-flankenpunkt_unten+gamma]),
										sphere_xlate([rg_outside, start*1.001, spiegelpunkt-flankenpunkt_unten]),								
										sphere_xlate([rg_outside, delta_f, flankenpunkt_unten]),
										sphere_xlate([rg_innen, delta_f, flankenpunkt_unten+gamma]),
										sphere_xlate([rg_innen, delta_f, spiegelpunkt-flankenpunkt_unten+gamma]),
										sphere_xlate([rg_outside, delta_f, spiegelpunkt-flankenpunkt_unten])								
									],
									faces = [[0,1,2],[0,2,3],[0,4,1],[1,4,5],[1,5,2],[2,5,6],[2,6,3],[3,6,7],[0,3,7],[0,7,4],[4,6,5],[4,7,6]],
									convexity =1
								);
							}
							// Zahn
							for (delta = [start:a_step:delta_a-a_step]){
								flankenpunkt_unten = sphere_calc(delta_b, delta);
								flankenpunkt_oben = sphere_calc(delta_b, delta+a_step);
								polyhedron(
									points = [
										sphere_xlate([rg_outside, delta, flankenpunkt_unten]),
										sphere_xlate([rg_innen, delta, flankenpunkt_unten+gamma]),
										sphere_xlate([rg_innen, delta, spiegelpunkt-flankenpunkt_unten+gamma]),
										sphere_xlate([rg_outside, delta, spiegelpunkt-flankenpunkt_unten]),								
										sphere_xlate([rg_outside, delta+a_step, flankenpunkt_oben]),
										sphere_xlate([rg_innen, delta+a_step, flankenpunkt_oben+gamma]),
										sphere_xlate([rg_innen, delta+a_step, spiegelpunkt-flankenpunkt_oben+gamma]),
										sphere_xlate([rg_outside, delta+a_step, spiegelpunkt-flankenpunkt_oben])									
									],
									faces = [[0,1,2],[0,2,3],[0,4,1],[1,4,5],[1,5,2],[2,5,6],[2,6,3],[3,6,7],[0,3,7],[0,7,4],[4,6,5],[4,7,6]],
									convexity =1
								);
							}
						}
					}
				}	
			}
		}
	}
}

/*  Pfeil-bevel_gear; uses the module "bevel_gear"
    modul = Height of the tooth head above the pitch circle
    teeth = Number of wheel teeth
    cone_angle, gear_width
    bore_hole = Center hole diameter
    contact_angle = contact_angle, Standard value = 20° according to DIN 867. Should not be greater than 45° sein.
    helix_angle = Schrägungswinkel, Standardwert = 0° */
module arrow_bevel_gear(modul, teeth, cone_angle, gear_width, bore_hole, contact_angle = 20, helix_angle=0){

	// Dimension calculations
	
	gear_width = gear_width / 2;
	
	d_outside = modul * teeth;								// Pitch cone diameter on the cone base surface,
																// corresponds to the chord in the spherical section
	r_outside = d_outside / 2;									// Part of sphere radius auf der Kegelgrundfläche 
	rg_outside = r_outside/sin(cone_angle);					// Large cone radius, corresponds to the length of the cone flank
	c = modul / 6;												// Header play
	df_outside = d_outside - (modul +c) * 2 * cos(cone_angle);
	rf_outside = df_outside / 2;
	delta_f = asin(rf_outside/rg_outside);
	height_f = rg_outside*cos(delta_f);							// Höhe des Kegels vom Fußkegel

	// Torsionswinkel gamma aus Schrägungswinkel
	gamma_g = 2*atan(gear_width*tan(helix_angle)/(2*rg_outside-gear_width));
	gamma = 2*asin(rg_outside/r_outside*sin(gamma_g/2));
	
	echo("Pitch cone diameter on the cone base surface = ", d_outside);
	
	// Größen für Komplementär-Kegelstumpf
	height_k = (rg_outside-gear_width)/cos(cone_angle);		// Höhe des Komplementärkegels für richtige Zahnlänge
	rk = (rg_outside-gear_width)/sin(cone_angle);			// Fußradius des Komplementärkegels
	rfk = rk*height_k*tan(delta_f)/(rk+height_k*tan(delta_f));	// Kopfradius des Zylinders für 
																// Komplementär-Kegelstumpf
	height_fk = rk*height_k/(height_k*tan(delta_f)+rk);			// height des Komplementär-Kegelstumpfs
	
	modul_innen = modul*(1-gear_width/rg_outside);

		union(){
		bevel_gear(modul, teeth, cone_angle, gear_width, bore_hole, contact_angle, helix_angle);		// untere Hälfte
		translate([0,0,height_f-height_fk])
			rotate(a=-gamma,v=[0,0,1])
				bevel_gear(modul_innen, teeth, cone_angle, gear_width, bore_hole, contact_angle, -helix_angle);	// obere Hälfte
	}
}

/*  Spiral-bevel_gear; uses the module "bevel_gear"
    modul = Height of the tooth head above the pitch circle
    teeth = Number of wheel teeth
    height = Höhe des Zahnrads
    bore_hole = Center hole diameter
    contact_angle = contact_angle, Standard value = 20° according to DIN 867. Should not be greater than 45° sein.
    helix_angle = Schrägungswinkel, Standardwert = 0° */
module spiralbevel_gear(modul, teeth, cone_angle, gear_width, bore_hole, contact_angle = 20, helix_angle=30){

	a_stepe = 16;

	// Dimension calculations
	
	b = gear_width / a_stepe;	
	d_outside = modul * teeth;								// Pitch cone diameter on the cone base surface,
																// corresponds to the chord in the spherical section
	r_outside = d_outside / 2;									// Part of sphere radius auf der Kegelgrundfläche 
	rg_outside = r_outside/sin(cone_angle);					// Large cone radius, corresponds to the length of the cone flank
	rg_mitte = rg_outside-gear_width/2;

	echo("Pitch cone diameter on the cone base surface = ", d_outside);

	a=tan(helix_angle)/rg_mitte;
	
	union(){
	for(i=[0:1:a_stepe-1]){
		r = rg_outside-i*b;
		helix_angle = a*r;
		modul_r = modul-b*i/rg_outside;
		translate([0,0,b*cos(cone_angle)*i])
			
			rotate(a=-helix_angle*i,v=[0,0,1])
				bevel_gear(modul_r, teeth, cone_angle, b, bore_hole, contact_angle, helix_angle);	// obere Hälfte
		}
	}
}

/*	cone_wheel_pair mit beliebigem axis_angle; uses the module "bevel_gear"
    modul = Höhe des Zahnkopfes über dem Teilkegel; Angabe für die Aussenseite des Kegels
    teeth_rad = Number of wheel teeth am Rad
    teeth_pinion = Number of wheel teeth am Ritzel
	axis_angle = Winkel zwischen den Achsen von Rad und Ritzel
    gear_width = width der Zähne von der Außenseite in Richtung Kegelspitze
    bore_hole_rad = Center hole diameter des Rads
    bore_hole_pinion = Center hole diameteren des Ritzels
    contact_angle = contact_angle, Standard value = 20° according to DIN 867. Should not be greater than 45° sein.
	helix_angle = Schrägungswinkel, Standardwert = 0°
	assembled = Komponenten zusammengebaut für Konstruktion oder auseinander zum 3D-Druck */
module cone_wheel_pair(modul, teeth_rad, teeth_pinion, axis_angle=90, gear_width, bore_hole_rad, bore_hole_pinion, contact_angle=20, helix_angle=0, assembled=true){
 
	// Dimension calculations
	r_rad = modul*teeth_rad/2;							// Partial cone radius of the wheel
	delta_rad = atan(sin(axis_angle)/(teeth_pinion/teeth_rad+cos(axis_angle)));	// Kegelwinkel des Rads
	delta_pinion = atan(sin(axis_angle)/(teeth_rad/teeth_pinion+cos(axis_angle)));// Kegelwingel des Ritzels
	rg = r_rad/sin(delta_rad);								// Radius der Großkugel
	c = modul / 6;											// Header play
	df_pinion = pi*rg*delta_pinion/90 - 2 * (modul + c);	// Fußkegeldurchmesser auf der Großkugel 
	rf_pinion = df_pinion / 2;								// Base cone radius on the large sphere
	delta_f_pinion = rf_pinion/(pi*rg) * 180;				// Kopfkegelwinkel
	rkf_pinion = rg*sin(delta_f_pinion);					// Radius des Kegelfußes
	height_f_pinion = rg*cos(delta_f_pinion);				// Höhe des Kegels vom Fußkegel
	
	echo("Kegelwinkel Rad = ", delta_rad);
	echo("Kegelwinkel Ritzel = ", delta_pinion);
 
	df_rad = pi*rg*delta_rad/90 - 2 * (modul + c);			// Fußkegeldurchmesser auf der Großkugel 
	rf_rad = df_rad / 2;									// Base cone radius on the large sphere
	delta_f_rad = rf_rad/(pi*rg) * 180;						// Kopfkegelwinkel
	rkf_rad = rg*sin(delta_f_rad);							// Radius des Kegelfußes
	height_f_rad = rg*cos(delta_f_rad);						// Höhe des Kegels vom Fußkegel

	echo("Höhe Rad = ", height_f_rad);
	echo("Höhe Ritzel = ", height_f_pinion);
	
	drehen = is_straight(teeth_pinion);
	
	// Drawing
	// Rad
	rotate([0,0,180*(1-spiel)/teeth_rad*drehen])
		bevel_gear(modul, teeth_rad, delta_rad, gear_width, bore_hole_rad, contact_angle, helix_angle);
	
	// Ritzel
	if (assembled)
		translate([-height_f_pinion*cos(90-axis_angle),0,height_f_rad-height_f_pinion*sin(90-axis_angle)])
			rotate([0,axis_angle,0])
				bevel_gear(modul, teeth_pinion, delta_pinion, gear_width, bore_hole_pinion, contact_angle, -helix_angle);
	else
		translate([rkf_pinion*2+modul+rkf_rad,0,0])
			bevel_gear(modul, teeth_pinion, delta_pinion, gear_width, bore_hole_pinion, contact_angle, -helix_angle);
 }

/*	Pfeil-cone_wheel_pair mit beliebigem axis_angle; uses the module "arrow_bevel_gear"
    modul = Höhe des Zahnkopfes über dem Teilkegel; Angabe für die Aussenseite des Kegels
    teeth_rad = Number of wheel teeth am Rad
    teeth_pinion = Number of wheel teeth am Ritzel
	axis_angle = Winkel zwischen den Achsen von Rad und Ritzel
    gear_width = width der Zähne von der Außenseite in Richtung Kegelspitze
    bore_hole_rad = Center hole diameter des Rads
    bore_hole_pinion = Center hole diameteren des Ritzels
    contact_angle = contact_angle, Standard value = 20° according to DIN 867. Should not be greater than 45° sein.
    helix_angle = Schrägungswinkel, Standardwert = 0°
	assembled = Komponenten zusammengebaut für Konstruktion oder auseinander zum 3D-Druck */
module locking_wheel_pair(modul, teeth_rad, teeth_pinion, axis_angle=90, gear_width, bore_hole_rad, bore_hole_pinion, contact_angle = 20, helix_angle=10, assembled=true){
 
	r_rad = modul*teeth_rad/2;							// Partial cone radius of the wheel
	delta_rad = atan(sin(axis_angle)/(teeth_pinion/teeth_rad+cos(axis_angle)));	// Kegelwinkel des Rads
	delta_pinion = atan(sin(axis_angle)/(teeth_rad/teeth_pinion+cos(axis_angle)));// Kegelwingel des Ritzels
	rg = r_rad/sin(delta_rad);								// Radius der Großkugel
	c = modul / 6;											// Header play
	df_pinion = pi*rg*delta_pinion/90 - 2 * (modul + c);	// Fußkegeldurchmesser auf der Großkugel 
	rf_pinion = df_pinion / 2;								// Base cone radius on the large sphere
	delta_f_pinion = rf_pinion/(pi*rg) * 180;				// Kopfkegelwinkel
	rkf_pinion = rg*sin(delta_f_pinion);					// Radius des Kegelfußes
	height_f_pinion = rg*cos(delta_f_pinion);				// Höhe des Kegels vom Fußkegel
	
	echo("Kegelwinkel Rad = ", delta_rad);
	echo("Kegelwinkel Ritzel = ", delta_pinion);
 
	df_rad = pi*rg*delta_rad/90 - 2 * (modul + c);			// Fußkegeldurchmesser auf der Großkugel 
	rf_rad = df_rad / 2;									// Base cone radius on the large sphere
	delta_f_rad = rf_rad/(pi*rg) * 180;						// Kopfkegelwinkel
	rkf_rad = rg*sin(delta_f_rad);							// Radius des Kegelfußes
	height_f_rad = rg*cos(delta_f_rad);						// Höhe des Kegels vom Fußkegel

	echo("Höhe Rad = ", height_f_rad);
	echo("Höhe Ritzel = ", height_f_pinion);
	
	drehen = is_straight(teeth_pinion);
	
	// Rad
	rotate([0,0,180*(1-spiel)/teeth_rad*drehen])
		arrow_bevel_gear(modul, teeth_rad, delta_rad, gear_width, bore_hole_rad, contact_angle, helix_angle);
	
	// Ritzel
	if (assembled)
		translate([-height_f_pinion*cos(90-axis_angle),0,height_f_rad-height_f_pinion*sin(90-axis_angle)])
			rotate([0,axis_angle,0])
				arrow_bevel_gear(modul, teeth_pinion, delta_pinion, gear_width, bore_hole_pinion, contact_angle, -helix_angle);
	else
		translate([rkf_pinion*2+modul+rkf_rad,0,0])
			arrow_bevel_gear(modul, teeth_pinion, delta_pinion, gear_width, bore_hole_pinion, contact_angle, -helix_angle);

}

/*
Berechnet eine screw / archimedische Schraube.
modul = Höhe des screwnkopfes über dem Teilzylinder
gear_number = Anzahl der Gänge (Zähne) der screw
length = Länge der screw
bore_hole = Center hole diameter
contact_angle = contact_angle, Standard value = 20° according to DIN 867. Should not be greater than 45° sein.
pitch_angle = pitch_angle der screw, entspricht 90° minus Schrägungswinkel. Positiver pitch_angle = rechtsdrehend.
assembled = Komponenten zusammengebaut für Konstruktion oder auseinander zum 3D-Druck */
module screw(modul, gear_number, length, bore_hole, contact_angle=20, pitch_angle, assembled=true){

	// Dimension calculations
	c = modul / 6;												// Header play
	r = modul*gear_number/(2*sin(pitch_angle));				// Teilzylinder-Radius
	rf = r - modul - c;											// Fußzylinder-Radius
	a = modul*gear_number/(90*tan(contact_angle));				// Spiralparameter
	tau_max = 180/gear_number*tan(contact_angle);				// Winkel von Fuß zu Kopf in der Normalen
	gamma = -rad*length/((rf+modul+c)*tan(pitch_angle));	// Torsionswinkel für Extrusion
	
	a_step = tau_max/16;
	
	// Drawing: extrudiere mit Verwindung eine Flaeche, die von zwei archimedischen Spiralen eingeschlossen wird
	if (assembled) {
		rotate([0,0,tau_max]){
			linear_extrude(height = length, center = false, convexity = 10, twist = gamma){
				difference(){
					union(){
						for(i=[0:1:gear_number-1]){
							polygon(
								concat(							
									[[0,0]],
									
									// ansteigende Zahnflanke
									[for (tau = [0:a_step:tau_max])
										pole2point([spirale(a, rf, tau), tau+i*(360/gear_number)])],
										
									// Zahnkopf
									[for (tau = [tau_max:a_step:180/gear_number])
										pole2point([spirale(a, rf, tau_max), tau+i*(360/gear_number)])],
									
									// absteigende Zahnflanke
									[for (tau = [180/gear_number:a_step:(180/gear_number+tau_max)])
										pole2point([spirale(a, rf, 180/gear_number+tau_max-tau), tau+i*(360/gear_number)])]
								)
							);
						}
						circle(rf);
					}
					circle(bore_hole/2); // Mittelbore_hole
				}
			}
		}
	}
	else {
		difference(){
			union(){
				translate([1,r*1.5,0]){
					rotate([90,0,90])
						screw(modul, gear_number, length, bore_hole, contact_angle, pitch_angle, assembled=true);
				}
				translate([length+1,-r*1.5,0]){
					rotate([90,0,-90])
						screw(modul, gear_number, length, bore_hole, contact_angle, pitch_angle, assembled=true);
					}
				}
			translate([length/2+1,0,-(r+modul+1)/2]){
					cube([length+2,3*r+2*(r+modul+1),r+modul+1], center = true);
				}
		}
	}
}

/*
Calculates a screw_gear_set. The screwn wheel is a common spur_gear without globoid geometry.
modul = Höhe des screwnkopfes über dem Teilzylinder bzw. des Zahnkopfes über dem Teilkreis
teeth = Number of wheel teeth
gear_number = Anzahl der Gänge (Zähne) der screw
width = gear_width
length = Länge der screw
bore_hole = Diameter of the center bore_hole of the screw
bore_hole_rad = Center hole diameter des spur_gears
contact_angle = contact_angle, Standard value = 20° according to DIN 867. Should not be greater than 45° sein.
pitch_angle = pitch_angle der screw, entspricht 90°-Schrägungswinkel. Positiver pitch_angle = rechtsdrehend.
optimized = Holes to save material/weight
assembled =  Komponenten zusammengebaut für Konstruktion oder auseinander zum 3D-Druck */
module screw_gear_set(modul, teeth, gear_number, width, length, bore_hole, bore_hole_rad, contact_angle=20, pitch_angle, optimized=true, assembled=true){
	
	c = modul / 6;												// Header play
	r_screw = modul*gear_number/(2*sin(pitch_angle));		// Teilzylinder-Radius screw
	r_rad = modul*teeth/2;									// Part of sphere radius spur_gear
	rf_screw = r_screw - modul - c;						// Fußzylinder-Radius
	gamma = -90*width*sin(pitch_angle)/(pi*r_rad);			// Rotationswinkel spur_gear
	tooth_space = modul*pi/cos(pitch_angle);				// tooth_space im Transversalschnitt
	x = is_straight(gear_number)? 0.5 : 1;

	if (assembled) {
		translate([r_screw,(ceil(length/(2*tooth_space))-x)*tooth_space,0])
			rotate([90,180/gear_number,0])
				screw(modul, gear_number, length, bore_hole, contact_angle, pitch_angle, assembled);

		translate([-r_rad,0,-width/2])
			rotate([0,0,gamma])
				spur_gear (modul, teeth, width, bore_hole_rad, contact_angle, -pitch_angle, optimized);
	}
	else {	
		screw(modul, gear_number, length, bore_hole, contact_angle, pitch_angle, assembled);

		translate([-2*r_rad,0,0])
			spur_gear (modul, teeth, width, bore_hole_rad, contact_angle, -pitch_angle, optimized);
	}
}

// linear_rack(modul=1, length=30, height=5, width=5, contact_angle=20, helix_angle=20);

// spur_gear (modul=1, teeth=30, width=5, bore_hole=4, contact_angle=20, helix_angle=20, optimized=true);

// arrow_gear (modul=1, teeth=30, width=5, bore_hole=4, contact_angle=20, helix_angle=30, optimized=true);

// rack_and_wheel (modul=1, length_stange=50, teeth_rad=30, height_stange=4, bore_hole_rad=4, width=5, contact_angle=20, helix_angle=0, assembled=true, optimized=true);

// angled_inner_gear (modul=1, teeth=30, width=5, randwidth=3, contact_angle=20, helix_angle=20);

// arrow_inner_gear (modul=1, teeth=30, width=5, randwidth=3, contact_angle=20, helix_angle=30);

// planetary_gear(modul=1, teeth_sonne=16, teeth_planet=9, anzahl_planeten=5, width=5, randwidth=3, bore_hole=4, contact_angle=20, helix_angle=30, assembled=true, optimized=true);

// bevel_gear(modul=1, teeth=30,  cone_angle=45, gear_width=5, bore_hole=4, contact_angle=20, helix_angle=20);

// arrow_bevel_gear(modul=1, teeth=30, cone_angle=45, gear_width=5, bore_hole=4, contact_angle=20, helix_angle=30);

// cone_wheel_pair(modul=1, teeth_rad=30, teeth_pinion=11, axis_angle=100, gear_width=5, bore_hole=4, contact_angle = 20, helix_angle=20, assembled=true);

// cone_wheel_pair(modul=1, teeth_rad=30, teeth_pinion=11, axis_angle=100, gear_width=5, bore_hole_rad=3, bore_hole_pinion=3, contact_angle=20, helix_angle=20, assembled=true);

// locking_wheel_pair(modul=1, teeth_rad=30, teeth_pinion=11, axis_angle=100, gear_width=5, bore_hole_rad=3, bore_hole_pinion=3, contact_angle = 20, helix_angle=30, assembled=false);

// screw(modul=1, gear_number=2, length=15, bore_hole=4, contact_angle=20, pitch_angle=10, assembled=true);

// screw_gear_set(modul=1, teeth=30, gear_number=2, width=8, length=20, bore_hole=4, bore_hole_rad=4, contact_angle=20, pitch_angle=10, optimized=true, assembled=true);