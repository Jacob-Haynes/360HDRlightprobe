INSTRUCTIONS

default calibration values: (can just run this section)
lower.ss = [ -9.343383e+02 0.000000e+00 2.546344e-04 -2.450943e-08 9.461811e-11  ];
lower.invpol = [ 1481.140014 812.598208 -58.068926 152.720328 50.274832 -4.348677 50.556259 -24.269045 -10.887808 52.821962 3.440259 -29.792634 -2.530701 8.810030 2.566491 ];
lower.xc = 1725.904819;
lower.yc = 1728.486293;
lower.c = 0.993223;
lower.d = 0.000077;
lower.e = -0.000119;
lower.height = 3456;
lower.width = 3456;
upper.ss = [-8.706589e+02 0.000000e+00 3.503911e-04 -1.414086e-07 1.222870e-10 ];
upper.invpol = [ 1464.982183 885.649072 -82.652189 119.590497 99.176419 -53.771453 58.137939 43.968097 -84.205343 7.232620 77.777963 -6.422891 -38.232782 -2.010324 9.544693 2.576489 ];
upper.xc = 1725.066534;
upper.yc = 1725.747413;
upper.c = 0.991498;
upper.d = 0.001621;
upper.e = 0.000991;
upper.height = 3456;
upper.width = 3456;
R = [-0.9912586, -0.0197387, 0.1304486; 0.0163975, -0.9995107, -0.0266381; 0.1309106, -0.0242662, 0.9910972];

To run the calibration for R matrix type: calibrate_R
follow comandline instructions

To form a 360 image from any raw (.tiff) file type: 
merge360(imported image, calibration lower, calibration upper, calibrated R)
