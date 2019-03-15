# Vect
c++ general purpose multithreaded numerical vector class
implements: operators, scalar vector aritmetics, several creators, memory mgr.
seq, reverse, apply func, filters, . etc.
uses kahan algorith to preserve precission when summing vector's components.
samples:

  VectReal vx(0, M_PI, 1e-4), vy(vx, sin), vytest=vx.func(sin);
  assert(vx.func(tan) == vx.func(sin) / vx.func(cos)); 
