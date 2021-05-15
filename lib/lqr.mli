open Owl

val classical_lqr : r2:float -> q:Mat.mat -> a:Mat.mat -> b:Mat.mat -> Mat.mat

val output_lqr
  :  r2:float
  -> q:Mat.mat
  -> a:Mat.mat
  -> b:Mat.mat
  -> param:(int * int) * Mat.mat
  -> Mat.mat
