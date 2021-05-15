open Owl

module type G = sig
  val m : Mat.mat
  val m_norm : Mat.mat
  val evals : Mat.mat
  val evecs : Mat.mat
  val get : int -> Mat.mat
  val top : int -> Mat.mat
  val bottom : int -> Mat.mat
end

module type T = sig
  val n : int

  module O : G
  module C : G
  module OSub : G
  module CSub : G
end

module type Prms = sig
  val a : Mat.mat
  val m1_slice : int list
  val gamma : Mat.mat option
end

module Make (P : Prms) : T
