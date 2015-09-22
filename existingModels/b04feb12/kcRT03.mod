: $Id: kcRT03.mod,v 1.5 2003/09/30 21:11:43 billl Exp $
TITLE Potasium C type current for RD Traub, J Neurophysiol 89:909-921, 2003

COMMENT

	Implemented by Maciej Lazarewicz 2003 (mlazarew@seas.upenn.edu)

ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
}
 
NEURON { 
	SUFFIX kcRT03
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE  gmax, ik, g, i
        GLOBAL alpha,beta
}

PARAMETER { 
  gmax = 0.0 	(mho/cm2)
  v ek 		(mV)  
  cai		(1)
} 

ASSIGNED { 
  g i
  ik 		(mA/cm2) 
  alpha beta	(/ms)
}
 
STATE {
  m
}

BREAKPOINT { 
  SOLVE states METHOD cnexp
  if( 0.004 * cai < 1 ) {
    g = gmax * m * 0.004 * cai
  } else {
    g = gmax * m
  }
  i = g * (v-ek) 
  ik=i
}
 
INITIAL { 
	settables(v) 
	m = alpha / ( alpha + beta )
	m = 0
}
 
DERIVATIVE states { 
	settables(v) 
	m' = alpha * ( 1 - m ) - beta * m 
}

UNITSOFF 

PROCEDURE settables(v) { 

	if( v < -10.0 ) {
		alpha = 2 / 37.95 * ( exp( ( v + 50 ) / 11 - ( v + 53.5 ) / 27 ) )

		: Note that there is typo in the paper - missing minus sign in the front of 'v'
		beta  = 2 * exp( ( - v - 53.5 ) / 27 ) - alpha
	}else{
		alpha = 2 * exp( ( - v - 53.5 ) / 27 )
		beta  = 0
	}
}

UNITSON
