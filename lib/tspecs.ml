module type T = sig
  val off1 : float
  val on : float
  val off2 : float
  val duration : float
  val gain_f : float -> float
  val input_gain_f : float -> float
  val sampling_dt : float
  val mov_input : float -> float
end

module Default : T = struct
  let off1 = Cmdargs.(get_float "-spont" |> default 0.4)
  let on = Cmdargs.(get_float "-prep" |> default 0.6)
  let off2 = Cmdargs.(get_float "-mov" |> default 0.6)
  let duration = off1 +. on +. off2
  let sampling_dt = 1E-3
  let gain_f = Signals.off_on_off ~off1 ~on ~off2
  let input_gain_f = gain_f
  let mov_input t = if t < off1 +. on then 0. else Defaults.mov_input (t -. off1 -. on)
end
