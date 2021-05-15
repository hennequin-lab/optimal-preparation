(* simple signals *)

let on_off_square ~on ~off =
  let period = on +. off in
  let rec f t = if t > period then f (t -. period) else if t < on then 1. else 0. in
  f


let off_on_square ~off ~on =
  let period = on +. off in
  let rec f t = if t > period then f (t -. period) else if t < off then 0. else 1. in
  f


let off_on_off ~off1 ~on ~off2 =
  let period = off1 +. on +. off2 in
  let rec f t =
    if t > period
    then f (t -. period)
    else if t < off1
    then 0.
    else if t -. off1 < on
    then 1.
    else 0.
  in
  f


let on_off_smooth ~tau_on ~tau_off ~on ~off =
  let period = on +. off in
  let max = 1. -. exp (-.on /. tau_on) in
  let rec f t =
    if t > period
    then f (t -. period)
    else if t < on
    then 1. -. exp (-.t /. tau_on)
    else max *. exp (-.(t -. on) /. tau_off)
  in
  f


let off_on_off_smooth ~tau_on ~tau_off ~off1 ~on ~off2 =
  let period = off1 +. on +. off2 in
  let max = 1. -. exp (-.on /. tau_on) in
  let rec f t =
    if t > period
    then f (t -. period)
    else if t < off1
    then 0.
    else if t -. off1 < on
    then 1. -. exp (-.(t -. off1) /. tau_on)
    else max *. exp (-.(t -. off1 -. on) /. tau_off)
  in
  f


let alpha_bump ~tau_rise ~tau_decay ~amplitude =
  let tmax =
    log (tau_decay /. tau_rise) *. tau_decay *. tau_rise /. (tau_decay -. tau_rise)
  in
  let f t = exp (-.t /. tau_decay) -. exp (-.t /. tau_rise) in
  let amax = f tmax in
  fun t -> amplitude /. amax *. f t
