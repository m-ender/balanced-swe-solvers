solver_unbalanced: rp1_shallow_roe_with_efix_unbalanced.f90
	f2py -m shallow_roe_with_efix_unbalanced -c $^

solver_leveque: rp1_shallow_roe_with_efix_leveque.f90
	f2py -m shallow_roe_with_efix_leveque -c $^

solver_rogers: rp1_shallow_roe_with_efix_rogers.f90
	f2py -m shallow_roe_with_efix_rogers -c $^

solver_rogers_geo: rp1_shallow_roe_with_efix_rogers_geo.f90
	f2py -m shallow_roe_with_efix_rogers_geo -c $^

all: solver_unbalanced solver_leveque solver_rogers
