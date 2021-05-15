open Owl
include Defaults_base

(* --------------------------------------------------------------------------------
   --     GOOD DEFAULT PARAMETERS                                                --
   -------------------------------------------------------------------------------- *)
module Main = struct
  include Base

  let w_rec = Mat.(load_txt (Dir.in_dir "w_rec"))
  let n = Mat.row_num w_rec
  let n_e = Maths.round (0.8 *. float n) |> int_of_float
  let n_i = n - n_e
  let a = Mat.(w_rec - eye n)
  let m1_slice = Array.init n_e (fun i -> i) |> Array.to_list
  let pmd_slice = Array.init n (fun i -> i) |> Array.to_list

  let gamma =
    let e = Mat.eye n in
    Mat.(e.!{L m1_slice; R []})


  let nl = Mat.relu
end
