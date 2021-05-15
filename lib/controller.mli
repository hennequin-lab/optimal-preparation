open Owl

(* All these controllers return thal_m1, m1_thal *)

(* here, thal_m1 = identity, m1_thal = optimal LQR feedback *)
val vanilla
  :  a:Mat.mat
  -> pmd_slice:int list
  -> r2:float
  -> q:Mat.mat
  -> Mat.mat * Mat.mat

(* dynamic feedback based on E neurons only 
   returns an N x N_E matrix *)
val dynamic
  :  a:Mat.mat
  -> pmd_slice:int list
  -> n_e:int
  -> tau:float
  -> taus:float * float
  -> r2:float
  -> q:Mat.mat
  -> Mat.mat * Mat.mat

(* decompose a loop into the product of 
   1. a realistic E/I B of size (n, n_z)
   2. a purely positive K of size (n_z, n_e)
   we regularize with the squared Frobenius norm of both matrices
   and eventually equalize their Frobenius norms after decomposition

   return accuracy, (K_xz, K_zy, K_yx) *)
val decompose
  :  n:int
  -> n_z:int
  -> pmd_slice:int list
  -> reg:float
  -> Mat.mat * Mat.mat
  -> float * (Mat.mat * Mat.mat * Mat.mat)
