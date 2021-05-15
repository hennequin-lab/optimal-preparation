open Owl
open Miscellaneous

(* --------------------------------------------------------------------------------
   --     GRAMIAN MATRICES, CONTROLLABLE/OBSERVABLE SUBSPACES                    --
   -------------------------------------------------------------------------------- *)

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

module Make (P : Prms) : T = struct
  open P

  let n, _ = Mat.shape a

  let get trans rhs =
    let rhs =
      match rhs with
      | `Identity -> Mat.(neg (eye n))
      | `Matrix ctc -> Mat.(neg ctc)
    in
    let z =
      Linalg.D.lyapunov
        (match trans with
        | `N -> a
        | `T -> Mat.transpose a)
        rhs
    in
    let z_norm = Mat.(float n /. trace z $* z) in
    let evecs, evals, _ = Linalg.D.svd z in
    let get i = Mat.col evecs i in
    let top k = Mat.cols evecs (Array.init k identity) in
    let bottom k = Mat.cols evecs (Array.init k (( + ) (n - k))) in
    z, z_norm, Mat.transpose evals, evecs, get, top, bottom


  module O = struct
    let m, m_norm, evals, evecs, get, top, bottom = get `T `Identity
  end

  module C = struct
    let m, m_norm, evals, evecs, get, top, bottom = get `N `Identity
  end

  module OSub = struct
    let m, m_norm, evals, evecs, get, top, bottom =
      let rhs =
        match gamma with
        | Some gamma -> Mat.(transpose gamma *@ gamma)
        | None -> Mat.(eye n)
      in
      get `T (`Matrix rhs)
  end

  module CSub = struct
    let m =
      let z = Linalg.D.lyapunov a Mat.(neg OSub.m_norm) in
      Mat.(z.!{L m1_slice; L m1_slice})


    let m_norm = Mat.(float (List.length m1_slice) /. trace m $* m)
    let evecs, evals, _ = Linalg.D.svd m
    let evals = Mat.transpose evals
    let get i = Mat.col evecs i
    let top k = Mat.cols evecs (Array.init k identity)
    let bottom k = Mat.cols evecs (Array.init k (( + ) (n - k)))
  end
end
