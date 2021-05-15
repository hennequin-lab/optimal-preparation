let angle_i = -36.
let angle_f = 180. +. 36.
let radius_i = 0.05
let radius_f = 0.2

let reach_angles ?(angle_i = angle_i) ?(angle_f = angle_f) n_angles =
  Owl.Mat.linspace angle_i angle_f n_angles |> Owl.Mat.to_array


let reach_radius ?(radius_i = radius_i) ?(radius_f = radius_f) n_radius =
  Owl.Mat.linspace radius_i radius_f n_radius |> Owl.Mat.to_array
