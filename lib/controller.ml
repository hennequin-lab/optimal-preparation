open Printf
open Owl

(* returns an N x N matrix *)
let vanilla ~a ~pmd_slice ~r2 ~q =
  let n, _ = Mat.shape a in
  let n_pmd = List.length pmd_slice in
  let b = Mat.zeros n n_pmd in
  Mat.(b.!{L pmd_slice; R []} <- Mat.eye n_pmd);
  Mat.eye n, Mat.(b *@ Lqr.classical_lqr ~r2 ~q ~a ~b)


let dynamic ~a ~pmd_slice ~n_e ~tau ~taus ~r2 ~q =
  let n, _ = Mat.shape a in
  let n_pmd = List.length pmd_slice in
  let n_i = n - n_e in
  let tau_z, tau_y = taus in
  let const_y = Mat.(tau /. tau_y $* eye n_e) in
  let const_z = Mat.(tau /. tau_z $* eye n_e) in
  let yx =
    (* sparsity p, and non-zero weight drawn from a lognormal distribution *)
    let p = 0.1 in
    let x = Mat.(bernoulli ~p n_e n_e /$ (p *. float n_e)) in
    Mat.(concat_horizontal x (zeros n_e n_i))
  in
  let a =
    Mat.(
      concat_vh
        [| [| a; zeros n n_e; zeros n n_e |]
         ; [| tau /. tau_y $* yx; neg const_y; zeros n_e n_e |]
         ; [| zeros n_e n; const_z; neg const_z |]
        |])
  in
  let b = Mat.zeros (n + n_e + n_e) n_pmd in
  Mat.(b.!{L pmd_slice; R []} <- Mat.eye n_pmd);
  let c = Mat.(concat_vh [| [| zeros n_e n; zeros n_e n_e; eye n_e |] |]) in
  let q =
    Mat.(
      concat_vh
        [| [| q; zeros n n_e; zeros n n_e |]
         ; [| zeros n_e n; zeros n_e n_e; zeros n_e n_e |]
         ; [| zeros n_e n; zeros n_e n_e; zeros n_e n_e |]
        |])
  in
  let xy = Lqr.output_lqr ~r2 ~q ~a ~b ~param:((n_pmd, n_e), c) in
  yx, Mat.(b *@ xy)


(* k_xy is further decomposed into k_xz * k_zy with appropriate
                       sign constraints, by the function [decompose] below *)

(* decompose a loop into the product of
   1. a realistic E/I B of size (n, n_z)
   2. a purely positive K of size (n_z, n_e)
   we regularize with the squared Frobenius norm of both matrices
   and eventually equalize their Frobenius norms after decomposition *)
let decompose ~n ~n_z ~pmd_slice ~reg (yx, xy) =
  let open Bigarray in
  let xy = Mat.(xy.!{L pmd_slice; R []}) in
  let n_pmd, n_e = Mat.shape xy in
  let n_z_e = int_of_float (Maths.round (0.5 *. float n_z)) in
  let n_prms = (n_pmd * n_z) + (n_z * n_e) in
  let lambda = 1. /. Mat.l2norm_sqr' xy in
  let reg = reg /. float n_prms in
  let for_zy prms =
    reshape_2 (genarray_of_array1 (Array1.sub prms 0 (n_z * n_e))) n_z n_e
  in
  let for_xz prms =
    reshape_2 (genarray_of_array1 (Array1.sub prms (n_z * n_e) (n_pmd * n_z))) n_pmd n_z
  in
  let f =
    let zy_prms = reshape_2 (Mat.empty n_z n_e) n_z n_e in
    let xz_prms = reshape_2 (Mat.empty n_pmd n_z) n_pmd n_z in
    fun prms gradient ->
      Array2.blit (for_zy prms) zy_prms;
      Array2.blit (for_xz prms) xz_prms;
      let open Algodiff.D in
      let t = tag () in
      let zy_prms = make_reverse (Arr (genarray_of_array2 zy_prms)) t in
      let xz_prms = make_reverse (Arr (genarray_of_array2 xz_prms)) t in
      let zy = Maths.sqr zy_prms in
      let xz =
        (* E/I parameterisation *)
        let e_part = Maths.get_slice [ []; [ 0; pred n_z_e ] ] xz_prms in
        let i_part = Maths.get_slice [ []; [ n_z_e; pred n_z ] ] xz_prms in
        Maths.(concat ~axis:1 (sqr e_part) (neg (sqr i_part)))
      in
      let cost =
        (* forward pass *)
        Maths.(
          (F lambda * l2norm_sqr' ((xz *@ zy) - Arr xy))
          + (F reg * (l2norm_sqr' zy + l2norm_sqr' xz)))
      in
      (* reverse pass if necessary *)
      (match gradient with
      | Some g ->
        reverse_prop (F 1.) cost;
        (* reverse-mode diff *)
        (* pack the gradients *)
        Bigarray.Genarray.blit
          (unpack_arr (adjval zy_prms))
          (genarray_of_array2 (for_zy g));
        Bigarray.Genarray.blit
          (unpack_arr (adjval xz_prms))
          (genarray_of_array2 (for_xz g))
      | None -> ());
      unpack_flt cost, unpack_arr (primal zy), unpack_arr (primal xz)
  in
  (* this is what L-BFGS understands *)
  let f_df x g =
    let cost, _, _ = f x (Some g) in
    cost
  in
  (* initial parameter vector *)
  let prms = Bigarray.reshape_1 (Mat.gaussian ~sigma:0.1 n_prms 1) n_prms in
  let every_iter_do k cost =
    if k mod 10 = 0
    then (
      Gc.major ();
      let ac, _, _ = f prms None in
      printf "\riteration %6i | cost = %10.5f | accuracy = %10.5f%!" k cost ac)
  in
  let stop st =
    let k = Lbfgs.iter st in
    let cost = Lbfgs.previous_f st in
    every_iter_do (k - 1) cost;
    false
  in
  Lbfgs.(C.min ~print:No ~pgtol:0. ~factr:1E9 ~corrections:20 ~stop f_df prms |> ignore);
  let accuracy, zy, xz = f prms None in
  (* finally equalize the Frob. norms *)
  let __zy = Mat.l2norm' zy
  and __xz = Mat.l2norm' xz
  and __yx = Mat.l2norm' yx in
  let zy = Mat.(Maths.(sqrt (__xz *. __yx) /. __zy) $* zy)
  and xz = Mat.(Maths.(sqrt (__zy *. __yx) /. __xz) $* xz)
  and yx = Mat.(Maths.(sqrt (__zy *. __xz) /. __yx) $* yx) in
  let full_xz = Mat.zeros n Mat.(col_num xz) in
  Mat.(full_xz.!{L pmd_slice; R []} <- xz);
  accuracy, (full_xz, zy, yx)
