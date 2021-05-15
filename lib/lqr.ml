open Printf
open Owl

(* ----------------------------------------------------------------------------
   ---                   VANILLA FULL-STATE LQR                             ---
   ---------------------------------------------------------------------------- *)

(* solve the LQR problem - cost function is infinite-horizon integral(xQx + uRu)dt *)
let classical_lqr ~r2 ~q ~a ~b =
  let _, m = Mat.shape b in
  let r = Mat.(r2 $* eye m) in
  let x = Linalg.D.care a b q r in
  Mat.(neg (transpose b *@ x /$ r2))


(* ----------------------------------------------------------------------------
   ---                   OUTPUT LQR                                         ---
   ---------------------------------------------------------------------------- *)

(* solve our output LQR problem using gradient descent (BFGS) *)
let output_lqr ~r2 ~q ~a ~b ~param =
  (* r2 : input energy penalty (R = diag(r2))
     q : quadratic cost matrix
     param: (Z_dim, C) → parameterization of the controller as ZC with C fixed *)
  let n, _ = Mat.shape a in
  let id = Mat.eye n in
  let (z_rows, z_cols), c = param in
  let n_prms = z_rows * z_cols in
  let bt = Mat.transpose b in
  let ct = Mat.transpose c in
  let k_of z = Mat.(z *@ c) in
  (* Lagrangian; L(Z,P,S) = trace(P) + trace(gS)
     where g(P,Z) = A_cl' P + P A_cl  + Q + K(Z) R K(Z)'
     stationary point conditions:
     dL/dP = 0  →  A_cl S + S A_cl' + I = 0 (Lyapunov equation 1)
     dL/dS = 0  →  g(P,Z) = 0 (Lyapunov equation for P, for fixed Z)
     dL/dZ = 0  →  2 (P + RK(Z)) S C' = 0

     key point is: when S and P solve their respective equations, g=0 anyways,
                   so the gradient of trace(P) w.r.t Z subject to g=0 is exactly the last equation
                   where P and S solve their respective Lyapunov equations *)
  let s_of acl = Linalg.D.lyapunov acl (Mat.neg id) in
  let p_of acl k =
    Linalg.D.lyapunov (Mat.transpose acl) Mat.(neg (q + (r2 $* transpose k *@ k)))
  in
  let grad_z k p s = Mat.(2. $* bt *@ (p + (r2 $* b *@ k)) *@ s *@ ct) in
  (* start from a standard lqr based solution *)
  let z =
    let k = classical_lqr ~r2 ~q ~a ~b in
    Mat.(k *@ transpose (Linalg.D.linsolve (c *@ transpose c) c))
  in
  (* define the function to be optimised via L-BFGS;
     returns the value of the cost function + puts the gradient in g *)
  let f_df z g =
    let z =
      Bigarray.(genarray_of_array2 (reshape_2 (genarray_of_array1 z) z_rows z_cols))
    in
    let g =
      Bigarray.(genarray_of_array2 (reshape_2 (genarray_of_array1 g) z_rows z_cols))
    in
    let k = k_of z in
    let acl = Mat.(a + (b *@ k)) in
    let s = s_of acl in
    let p = p_of acl k in
    let g_ = grad_z k p s in
    Mat.copy_ g_ ~out:g;
    Mat.trace p
  in
  let stop iter cost =
    if iter mod 1 = 0
    then (
      let k = k_of z in
      let acl = Mat.(a + (b *@ k)) in
      let sa = Miscellaneous.spectral_abscissa acl in
      printf "\r iter %5i | SA = %.4f | trP = %.5f%!" iter sa cost);
    false
  in
  Lbfgs.(
    C.min
      ~print:No
      ~pgtol:1E-5
      ~factr:1E7
      ~corrections:10
      ~stop:(fun st -> stop (Lbfgs.iter st) (Lbfgs.previous_f st))
      f_df
      Bigarray.(reshape_1 z n_prms)
    |> ignore);
  print_newline ();
  z
