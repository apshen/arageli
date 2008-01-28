ARAGELICXXFLAGS = -DARAGELI_SOURCE -O2 -DARAGELI_DEBUG_LEVEL=1
TESTCXXFLAGS = -DARAGELI_SOURCE -O2 -g -DARAGELI_DEBUG_LEVEL=3 -DARAGELI_ASSERT_THROW_EXCEPTION

# Precompiled header(s)
arageli.hpp.gch: arageli.hpp
	$(CC) -x c++-header $(ARAGELICXXFLAGS) $< -o $@

# We are building some libraries. These libraries will be intalled in libdir.
lib_LIBRARIES = libarageli.a
libarageli_a_CXXFLAGS = $(ARAGELICXXFLAGS)
libarageli_a_SOURCES = \
						advsimplmeth.cpp\
						algebraic.cpp\
						algebrslt.cpp\
						bareiss.cpp\
						berlekamp_hensel.cpp\
						bigar.cpp\
						big_const.cpp\
						big_float.cpp\
						big_int.cpp\
						cantor_zassenhaus.cpp\
						cone.cpp\
						exception.cpp\
						gauss.cpp\
						gcd.cpp\
						generic.cpp\
						gram_schmidt.cpp\
						groebner_basis.cpp\
						hermite.cpp\
						intalg.cpp\
						intconvex.cpp\
						intcount_barvinok.cpp\
						io.cpp\
						karatsuba.cpp\
						linopt.cpp\
						lll.cpp\
						logarithm.cpp\
						matrix.cpp\
						motzkin_burger.cpp\
						orthogonal.cpp\
						pollard_mult.cpp\
						polyhedron.cpp\
						polyhedron_relation.cpp\
						polynom.cpp\
						powerest.cpp\
						prime.cpp\
						rational.cpp\
						resultant.cpp\
						sideset.cpp\
						smith.cpp\
						smithpoly.cpp\
						solve_linsys.cpp\
						sparse_multipolynom.cpp\
						sparse_polynom.cpp\
						sturm.cpp\
						tex_cartesian2d.cpp\
						texout.cpp\
						timer.cpp\
						triangulation.cpp\
						type_opers.cpp\
						type_traits.cpp\
						_utility.cpp\
						vecalg.cpp\
						vector.cpp\
						winograd.cpp\
						advsimplmeth.hpp\
						algebraic.hpp\
						algebrslt.hpp\
						arageli.hpp\
						auto_link.hpp\
						bareiss.hpp\
						basefuncs.hpp\
						berlekamp_hensel.hpp\
						bigar.hpp\
						big_const.hpp\
						big_float.hpp\
						big_int.hpp\
						cantor_zassenhaus.hpp\
						cmp.hpp\
						cone.hpp\
						config.hpp\
						ctrl_latexlog.hpp\
						ctrl_slog.hpp\
						exception.hpp\
						factory.hpp\
						fft.hpp\
						frwrddecl.hpp\
						functional.hpp\
						function_traits.hpp\
						gauss.hpp\
						gcd.hpp\
						generic.hpp\
						gram_schmidt.hpp\
						groebner_basis.hpp\
						hermite.hpp\
						intalg.hpp\
						intconvex.hpp\
						intcount_barvinok.hpp\
						interval.hpp\
						io.hpp\
						iomanip.hpp\
						iteradapt.hpp\
						karatsuba.hpp\
						linopt.hpp\
						lll.hpp\
						logarithm.hpp\
						matrix.hpp\
						misc.hpp\
						mixcomp.hpp\
						motzkin_burger.hpp\
						orthogonal.hpp\
						pollard_mult.hpp\
						polyalg.hpp\
						polyhedron.hpp\
						polyhedron_relation.hpp\
						polynom.hpp\
						powerest.hpp\
						prime.hpp\
						rand.hpp\
						random.hpp\
						rational.hpp\
						refcntr.hpp\
						residue.hpp\
						resultant.hpp\
						rndbase.hpp\
						setenumrnd.hpp\
						sideset.hpp\
						simplex_method.hpp\
						skeleton.hpp\
						smith.hpp\
						smithpoly.hpp\
						solve_linsys.hpp\
						sparse_multipolynom.hpp\
						sparse_polynom.hpp\
						sqrfree.hpp\
						std_import.hpp\
						sturm.hpp\
						tex_cartesian2d.hpp\
						texout.hpp\
						timer.hpp\
						triangulation.hpp\
						type_opers.hpp\
						type_pair_traits.hpp\
						type_traits.hpp\
						_utility.hpp\
						vecalg.hpp\
						vector.hpp\
						winograd.hpp

check_LIBRARIES = libtestarageli.a
libtestarageli_a_CXXFLAGS = $(TESTCXXFLAGS)
libtestarageli_a_SOURCES = $(libarageli_a_SOURCES)
