# Timings for the test suite
Taken from [Nasser's website](https://www.12000.org/my_notes/CAS_integration_tests/reports/rubi_4_16_1_graded/inch1.htm#x2-10001)

The current number of problems in this test suite is [71994].
1.1 Listing of CAS systems tested

The following systems were tested at this time.

    Mathematica 12.1 (64 bit) on windows 10.
    Rubi 4.16.1 in Mathematica 12 on windows 10.
    Maple 2020 (64 bit) on windows 10.
    Maxima 5.43 on Linux. (via sagemath 8.9)
    Fricas 1.3.6 on Linux (via sagemath 9.0)
    Sympy 1.5 under Python 3.7.3 using Anaconda distribution.
    Giac/Xcas 1.5 on Linux. (via sagemath 8.9)

System 	          solved 	       Failed
Rubi 	% 99.52 ( 71651 ) 	% 0.48 ( 343 )
Mathem 	% 98.35 ( 70804 ) 	% 1.65 ( 1190 )
Maple 	% 83.46 ( 60084 ) 	% 16.54 ( 11910 )
Fricas 	% 68.07 ( 49005 ) 	% 31.93 ( 22989 )
Giac 	% 52.51 ( 37804 ) 	% 47.49 ( 34190 )
Maxima 	% 43.03 ( 30982 ) 	% 56.97 ( 41012 )
Sympy 	% 32.41 ( 23332 ) 	% 67.59 ( 48662 )

System 	% A 	% B     % C     % F 
Rubi 	98.89 	0.23 	0.41 	0.48
Mathem  74.67 	6.18 	17.49 	1.65
Maple 	52.8 	22.93 	7.72 	16.54
Maxima 	33.24 	8.83 	0.96 	56.97
Fricas 	48.52 	18.03 	1.51 	31.93
Sympy 	25.08 	4.7 	2.63 	67.59
Giac 	39. 	12.5 	1.01 	47.49 

System 	Avg(s)  avgsize Norm.mean Med.size 	Norm.median
Rubi 	0.28 	156.75 	1. 	    107. 	1.
Mathem 	1.77 	800.29 	2.8 	92. 	0.94
Maple 	0.46 	62669. 	743.6 	131. 	1.27
Maxima 	1.34 	284.97 	2.46 	96. 	1.36
Fricas 	2.81 	935.28 	6.76 	302. 	3.43
Sympy 	9.65 	230.82 	2.53 	70. 	1.14
Giac 	1.53 	301.08 	2.55 	120. 	1.49 

72000 tests
Rubi 5.6 horus	
Mathem 35.4	
Maple LOL	
Maxima 	
Fricas 	
Sympy 	
Giac 	

Integrand type 	   probs 	Rubi 	Mathe 	Maple 	Maxima 	Fricas 	Sympy 	Giac
Independent tests 	1892 	98.31 	98.73 	92.18 	79.39 	94.34 	71.78 	82.72
Algebraic Binomial 	14276 	99.99 	99.7 	82.18 	42.02 	70.91 	59.27 	62.61
Algebraic Trinomia 	10187 	99.99 	98.89 	90.67 	38.56 	75.76 	40.39 	61.23
Algebraic Miscella 	1519 	98.62 	98.16 	87.23 	42.92 	74.06 	45.69 	54.18
Exponentials    	965 	99.17 	96.68 	80.21 	60.93 	87.67 	40.83 	46.74
Logarithms 	        3085 	98.51 	97.8 	54.49 	48.36 	57.76 	25.32 	43.37
Trigonometric 	    22551 	99.56 	97.61 	85.75 	41.34 	63.42 	13.64 	44.13
Inverse Trigonome 	4585 	99.65 	97.97 	83.84 	31.15 	48.29 	28.16 	48.05
Hyperbolic 	        5166 	98.32 	98.03 	82.58 	57.08 	84.84 	20.75 	62.45
Inverse Hyperbolic 	6626 	99.52 	98.46 	80.47 	40.34 	62.27 	24.89 	39.5
Special functions 	999 	100. 	95.4 	69.97 	35.54 	48.85 	39.34 	34.93 

